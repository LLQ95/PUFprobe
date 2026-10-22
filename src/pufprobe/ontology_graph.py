"""Gene Ontology graph and true-path (ontology) propagation.

The search channels report the most specific GO terms they can support. Gene
Ontology annotations, however, obey the true-path rule: a protein annotated to
a term is also annotated to every ancestor of that term along the ontology's
transitive ``is_a`` and ``part_of`` links. Without this closure a profile hit
to a specific term and a database hit to one of its parents cannot converge on
the same candidate, and CAFA-style evaluation requires predictions across the
whole ontology rather than only at the leaf terms.

This module parses the OBO 1.2 flat file shipped by the Gene Ontology
consortium (``go-basic.obo``; optionally gzip-compressed). ``go-basic.obo`` is
already filtered to the basic relations and guaranteed to be valid for the
true-path rule, so only ``is_a`` and ``part_of`` edges are followed upwards.
``regulates`` is deliberately not propagated: a process that regulates another
process is not an instance of it.

``alt_id`` terms are folded onto their primary id and obsolete terms are
recorded but excluded from the graph. The graph is intentionally small and
dependency-free so the rank stage can load it without a database server.
"""

from __future__ import annotations

import gzip
import re
from functools import lru_cache
from pathlib import Path
from typing import Dict, IO, List, Optional, Set, Tuple

from pufprobe.logging_utils import get_logger
from pufprobe.models import Evidence

logger = get_logger(__name__)

# Relations that imply the true-path rule when traversed toward the root.
# "is_a" is a first-class OBO tag; "part_of" arrives as a typed relationship.
_TRUE_PATH_RELATIONS = frozenset({"part_of"})

_IS_A_RE = re.compile(r"is_a:\s+(GO:\d{7})")
_REL_RE = re.compile(r"relationship:\s+(\w+)\s+(GO:\d{7})")
_ID_RE = re.compile(r"^id:\s+(GO:\d{7})\s*$")
_ALT_RE = re.compile(r"^alt_id:\s+(GO:\d{7})\s*$")
_NAME_RE = re.compile(r"^name:\s+(.+?)\s*$")
_NS_RE = re.compile(r"^namespace:\s+(\w+)\s*$")
_OBSOLETE_RE = re.compile(r"^is_obsolete:\s*true\s*$")


class GoOntology:
    """A directed acyclic GO graph keyed by primary GO identifiers."""

    def __init__(
        self,
        parents: Dict[str, Set[str]],
        names: Dict[str, str],
        namespaces: Dict[str, str],
        aliases: Dict[str, str],
        obsolete: Set[str],
    ) -> None:
        self._parents = parents
        self.names = names
        self.namespaces = namespaces
        self._aliases = aliases
        self.obsolete = obsolete

    @property
    def n_terms(self) -> int:
        return len(self._parents)

    def resolve(self, term_id: str) -> str:
        """Map an alt_id onto its primary id; pass through unknown ids."""
        return self._aliases.get(term_id, term_id)

    def name(self, term_id: str) -> str:
        """Human-readable label for a term, falling back to the identifier."""
        term_id = self.resolve(term_id)
        return self.names.get(term_id, term_id)

    def is_known(self, term_id: str) -> bool:
        term_id = self.resolve(term_id)
        return term_id in self._parents and term_id not in self.obsolete

    @lru_cache(maxsize=None)
    def ancestors(self, term_id: str) -> Tuple[str, ...]:
        """Return the term and all of its true-path ancestors.

        The term itself is included so callers can treat the returned set as
        the full closure of terms implied by an annotation. Unknown terms
        resolve to a closure containing only themselves.
        """
        term_id = self.resolve(term_id)
        if term_id in self.obsolete:
            return tuple()
        if term_id not in self._parents:
            # An unknown id has no parents to propagate to; keep it as an
            # isolated closure so the direct annotation is preserved.
            return (term_id,) if term_id.startswith("GO:") else tuple()
        closure: Set[str] = set()
        stack: List[str] = [term_id]
        while stack:
            current = stack.pop()
            if current in closure or current in self.obsolete:
                continue
            closure.add(current)
            stack.extend(self._parents.get(current, ()))
        return tuple(sorted(closure))

    def propagate_evidence(
        self, evidence_map: Dict[str, List[Evidence]]
    ) -> Tuple[Dict[str, List[Evidence]], int]:
        """Extend GO evidence with ancestor terms under the true-path rule.

        Non-GO evidence is returned unchanged. For every GO evidence object the
        full ancestor closure of its target is added in the same channel, so the
        propagated evidence keeps the same evidence group as the direct hit and
        never masquerades as an independent witness. For a given
        ``(ancestor, channel)`` only the strongest supporting evidence is kept;
        a direct annotation to the ancestor already present in the input is
        preferred over a propagated one when it is at least as strong.

        Returns the (new) evidence map and the number of added objects.
        """
        added = 0
        out_map: Dict[str, List[Evidence]] = {}

        for key, hits in evidence_map.items():
            go_hits = [h for h in hits if _is_go(h)]
            if not go_hits:
                out_map[key] = list(hits)
                continue

            # Strongest direct evidence already present for (term, channel);
            # these ancestors need no propagated duplicate when they are at
            # least as strong as the propagated support.
            direct_score: Dict[Tuple[str, str], float] = {}
            for hit in go_hits:
                term = self.resolve(hit.target_id)
                slot = (term, hit.channel.value)
                direct_score[slot] = max(direct_score.get(slot, float("-inf")), hit.score)

            propagated: Dict[Tuple[str, str], Evidence] = {}
            for hit in go_hits:
                term = self.resolve(hit.target_id)
                closure = self.ancestors(term)
                for ancestor in closure:
                    if ancestor == term or ancestor in self.obsolete:
                        continue
                    slot = (ancestor, hit.channel.value)
                    if direct_score.get(slot, float("-inf")) >= hit.score:
                        continue
                    current = propagated.get(slot)
                    if current is not None and current.score >= hit.score:
                        continue
                    propagated[slot] = self._propagated_evidence(hit, term, ancestor)

            out_map[key] = list(hits) + list(propagated.values())
            added += len(propagated)

        return out_map, added

    def _propagated_evidence(
        self, source: Evidence, source_term: str, ancestor: str
    ) -> Evidence:
        details = dict(source.details or {})
        details["db"] = "go"
        details["propagated"] = True
        details["propagated_from"] = source_term
        return Evidence(
            channel=source.channel,
            target_id=ancestor,
            target_label=self.name(ancestor),
            score=source.score,
            evalue=source.evalue,
            coverage=source.coverage,
            details=details,
        )


def _is_go(hit: Evidence) -> bool:
    tid = (hit.target_id or "").strip()
    return tid.startswith("GO:") or str((hit.details or {}).get("db", "")).lower() == "go"


def _open_obo(path: Path) -> IO[str]:
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def parse_go_obo(handle: IO[str]) -> GoOntology:
    """Parse an OBO 1.2 stream into a :class:`GoOntology`."""
    parents: Dict[str, Set[str]] = {}
    names: Dict[str, str] = {}
    namespaces: Dict[str, str] = {}
    aliases: Dict[str, str] = {}
    obsolete: Set[str] = set()

    current: Optional[Dict[str, object]] = None
    in_term = False

    def flush() -> None:
        nonlocal current
        if current is None:
            return
        term_id = current["id"]
        if current.get("obsolete"):
            obsolete.add(term_id)
            current = None
            return
        parents.setdefault(term_id, set()).update(current["parents"])  # type: ignore[arg-type]
        if current.get("name"):
            names[term_id] = str(current["name"])
        if current.get("namespace"):
            namespaces[term_id] = str(current["namespace"])
        for alt in current["alts"]:  # type: ignore[union-attr]
            aliases[alt] = term_id
        current = None

    for raw_line in handle:
        line = raw_line.rstrip("\n")
        if line.startswith("[Term]"):
            flush()
            in_term = True
            current = {"id": None, "parents": set(), "alts": [], "name": None,
                       "namespace": None, "obsolete": False}
            continue
        if line.startswith("["):
            # [Typedef], [Instance], ...: close any open term and ignore.
            flush()
            in_term = False
            continue
        if not in_term or current is None:
            continue

        match = _ID_RE.match(line)
        if match and current["id"] is None:
            current["id"] = match.group(1)
            continue
        match = _ALT_RE.match(line)
        if match:
            current["alts"].append(match.group(1))  # type: ignore[union-attr]
            continue
        match = _NAME_RE.match(line)
        if match and current["name"] is None:
            current["name"] = match.group(1)
            continue
        match = _NS_RE.match(line)
        if match and current["namespace"] is None:
            current["namespace"] = match.group(1)
            continue
        if _OBSOLETE_RE.match(line):
            current["obsolete"] = True
            continue
        match = _IS_A_RE.match(line)
        if match:
            current["parents"].add(match.group(1))  # type: ignore[union-attr]
            continue
        match = _REL_RE.match(line)
        if match and match.group(1) in _TRUE_PATH_RELATIONS:
            current["parents"].add(match.group(2))  # type: ignore[union-attr]

    flush()

    # Aliases must not remain dangling when their primary id lacked edges.
    for term_id in list(parents):
        parents.setdefault(term_id, set())
    return GoOntology(parents, names, namespaces, aliases, obsolete)


def load_go_obo(path: Optional[Path]) -> Optional[GoOntology]:
    """Load the GO graph from ``path``; return None when unavailable.

    A missing or unreadable file is logged and degrades to ``None`` so the
    rank stage keeps working without ontology propagation.
    """
    if path is None:
        return None
    path = Path(path)
    if not path.exists():
        logger.info("GO ontology not found at %s; true-path propagation disabled.", path)
        return None
    try:
        with _open_obo(path) as handle:
            ontology = parse_go_obo(handle)
    except (OSError, ValueError) as exc:
        logger.warning("Could not parse GO ontology %s (%s); propagation disabled.", path, exc)
        return None
    logger.info(
        "GO ontology loaded: %d terms, %d alt ids, %d obsolete (%s).",
        ontology.n_terms, len(ontology._aliases), len(ontology.obsolete), path,
    )
    return ontology
