#!/usr/bin/env python3
"""Build the cross-species hypothetical-protein family network.

Families are nodes and k-nearest-neighbour edges from the representative
embeddings connect structurally or sequence-similar families. Leiden
community detection groups the families into a handful of macro-communities,
which are the "clusters" reported in the paper. Each family and community is
annotated with size, taxonomic breadth, mean length and the fraction of
representatives that carry a structural/Pfam hit, so the largest and the best
annotated communities can be read directly from the tables.

Two resolutions are written for plotting: the full family-level graph (all
nodes with coordinates) and a community-level aggregate graph that shows the
macro-clusters and their relative sizes at a glance.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np


def load_embeddings(path: Path) -> Dict[str, np.ndarray]:
    vectors: Dict[str, np.ndarray] = {}
    if not path.exists():
        return vectors
    with np.load(path) as data:
        for key in data.files:
            vectors[key] = np.asarray(data[key], dtype=np.float32)
    return vectors


def load_families(path: Path) -> Dict[str, dict]:
    families: Dict[str, dict] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            families[row["family_id"]] = row
    return families


def load_members(path: Path) -> Dict[str, Counter]:
    """family_id -> Counter of phyla."""
    phyla: Dict[str, Counter] = defaultdict(Counter)
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            phyla[row["family_id"]][row["phylum"]] += 1
    return phyla


def load_structure_hits(path: Optional[Path]) -> Dict[str, dict]:
    hits: Dict[str, dict] = {}
    if path is None or not path.exists():
        return hits
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            hits[row["accession"]] = row
    return hits


def _rescale(coords) -> list:
    """Normalise coordinates to [-1, 1]; a degenerate axis maps to 0."""
    pts = [(float(x), float(y)) for x, y in coords]
    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]

    def norm(values):
        lo, hi = min(values), max(values)
        if hi - lo < 1e-6:
            return [0.0 for _ in values]
        return [(v - lo) / (hi - lo) * 2.0 - 1.0 for v in values]

    return list(zip(norm(xs), norm(ys)))


def build_knn(vectors: Dict[str, np.ndarray], keys: List[str], k: int
              ) -> List[tuple]:
    """Undirected kNN edges with cosine weights; duplicates keep the max."""
    from sklearn.neighbors import NearestNeighbors

    matrix = np.vstack([vectors[key] for key in keys])
    norms = np.linalg.norm(matrix, axis=1, keepdims=True)
    matrix = matrix / np.clip(norms, 1e-8, None)
    k_eff = min(k + 1, len(keys))
    nn = NearestNeighbors(n_neighbors=k_eff, metric="cosine", n_jobs=-1)
    nn.fit(matrix)
    distances, indices = nn.kneighbors(matrix)
    edge_map: Dict[tuple, float] = {}
    for i in range(len(keys)):
        for j, dist in zip(indices[i], distances[i]):
            if j == i:
                continue
            a, b = sorted((keys[i], keys[j]))
            weight = 1.0 - float(dist)
            if (a, b) not in edge_map or weight > edge_map[(a, b)]:
                edge_map[(a, b)] = weight
    return [(a, b, w) for (a, b), w in edge_map.items()]


def dominant_phylum(counter: Counter) -> tuple:
    if not counter:
        return "", 0.0
    phylum, count = counter.most_common(1)[0]
    return phylum, count / sum(counter.values())


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cluster-dir", type=Path, required=True)
    parser.add_argument("--structure-dir", type=Path, default=None)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--k", type=int, default=12)
    parser.add_argument("--resolution", type=float, default=1.0)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--prostt5", action="store_true",
                        help="concatenate ProstT5 embeddings with ESM-2")
    parser.add_argument("--core-min-size", type=int, default=2,
                        help="families at least this large are drawn in the detail plot")
    args = parser.parse_args()

    args.output.mkdir(parents=True, exist_ok=True)
    families = load_families(args.cluster_dir / "families.tsv")
    members_phyla = load_members(args.cluster_dir / "family_members.tsv")
    struct_hits = load_structure_hits(
        args.structure_dir / "structure_hits.tsv" if args.structure_dir else None)

    esm = load_embeddings(args.cluster_dir / "representative_embeddings.npz")
    if args.prostt5:
        t5 = load_embeddings(args.cluster_dir / "prostt5_embeddings.npz")
        if t5:
            shared = sorted(set(esm) & set(t5))
            if shared:
                esm = {key: np.concatenate([esm[key], t5[key]]) for key in shared}
                print(f"Combined ESM-2 + ProstT5 embeddings for {len(shared)} families.")

    # Map representative accession -> family id and keep families with vectors.
    rep_to_family = {fam["representative"]: fid for fid, fam in families.items()}
    keys = [acc for acc in sorted(esm) if acc in rep_to_family]
    if not keys:
        print("No representative embeddings available; run stage 08 first.",
              file=sys.stderr)
        return 1
    family_keys = [rep_to_family[acc] for acc in keys]
    print(f"Building kNN graph for {len(family_keys)} families (k={args.k}) ...")

    edges = build_knn(esm, keys, args.k)
    # convert accession endpoints to family ids for graph construction
    rep_edge = {(a, b): w for a, b, w in edges}
    family_edges = []
    for (a, b), w in rep_edge.items():
        family_edges.append((rep_to_family[a], rep_to_family[b], w))

    # --- igraph + Leiden -----------------------------------------------------
    import igraph as ig
    import leidenalg

    graph = ig.Graph.TupleList(
        [(a, b, w) for a, b, w in family_edges],
        directed=False, weights=True,
    )
    partition = leidenalg.find_partition(
        graph, leidenalg.RBConfigurationVertexPartition,
        weights="weight", n_iterations=-1,
        resolution_parameter=args.resolution, seed=args.seed,
    )
    community_of = {graph.vs[i]["name"]: part
                    for i, part in enumerate(partition.membership)}
    print(f"Leiden resolution {args.resolution} -> {len(partition)} communities.")

    # --- family-level node attributes ---------------------------------------
    def fam_size(fid: str) -> int:
        try:
            return int(families[fid]["size"])
        except (KeyError, ValueError):
            return 1

    family_phyla = {fid: members_phyla.get(fid, Counter()) for fid in family_keys}
    struct_pfam = {fid: (struct_hits.get(rep, {}).get("pfams", "") if rep in
                         struct_hits else "")
                   for fid, rep in ((fid, families[fid]["representative"])
                                    for fid in family_keys)}

    # Layout: force-directed on the core graph for the detail plot. Components
    # are laid out separately so a disconnected graph does not collapse onto
    # the origin.
    core_set = {fid for fid in family_keys if fam_size(fid) >= args.core_min_size}
    if not core_set:
        core_set = set(family_keys)
    layout_graph = graph.subgraph([v for v in graph.vs["name"] if v in core_set])
    core_xy = {}
    if layout_graph.vcount() > 0:
        try:
            base_layout = layout_graph.layout_fruchterman_reingold(weights="weight")
            coords = layout_graph.layout_components(base_layout)
        except Exception:  # noqa: BLE001 - fall back to a deterministic layout
            coords = layout_graph.layout_circle()
        coords = _rescale(coords)
        core_xy = {layout_graph.vs[i]["name"]: (coords[i][0], coords[i][1])
                   for i in range(layout_graph.vcount())}

    nodes_path = args.output / "nodes.tsv"
    with nodes_path.open("w", encoding="utf-8", newline="") as handle:
        header = ["family_id", "representative", "community", "size", "n_phyla",
                  "mean_length", "dominant_phylum", "phylum_fraction",
                  "structure_pfam", "core", "x", "y"]
        handle.write("\t".join(header) + "\n")
        for fid in family_keys:
            fam = families[fid]
            dom, frac = dominant_phylum(family_phyla.get(fid, Counter()))
            x, y = core_xy.get(fid, ("", ""))
            handle.write(
                f"{fid}\t{fam['representative']}\t{community_of.get(fid, '')}\t"
                f"{fam['size']}\t{fam['n_phyla']}\t{fam['mean_length']}\t"
                f"{dom}\t{frac:.3f}\t{struct_pfam.get(fid, '')}\t"
                f"{int(fid in core_set)}\t{x}\t{y}\n"
            )

    # --- community aggregate -------------------------------------------------
    members_of: Dict[int, List[str]] = defaultdict(list)
    for fid in family_keys:
        members_of[community_of.get(fid, 0)].append(fid)

    community_rows = []
    for part, fids in sorted(members_of.items(), key=lambda kv: -sum(
            fam_size(f) for f in kv[1])):
        n_families = len(fids)
        n_proteins = sum(fam_size(f) for f in fids)
        phyla_counter: Counter = Counter()
        for fid in fids:
            phyla_counter.update(family_phyla.get(fid, Counter()))
        dom, frac = dominant_phylum(phyla_counter)
        lengths = [float(families[f]["mean_length"]) for f in fids]
        cross_phylum = sum(1 for f in fids if int(families[f]["n_phyla"]) > 1)
        with_struct = sum(1 for f in fids if struct_pfam.get(f))
        median_length = float(np.median(lengths))
        label = _community_label(dom, frac, median_length, cross_phylum,
                                 n_families, with_struct)
        community_rows.append({
            "community": part,
            "n_families": n_families,
            "n_proteins": n_proteins,
            "dominant_phylum": dom,
            "phylum_fraction": round(frac, 3),
            "median_length": round(median_length, 1),
            "cross_phylum_families": cross_phylum,
            "structure_annotated_families": with_struct,
            "structure_annotated_fraction": round(with_struct / n_families, 3),
            "label": label,
        })

    # inter-community edges (aggregate weight)
    comm_edges: Counter = Counter()
    for a, b, w in family_edges:
        ca, cb = community_of.get(a), community_of.get(b)
        if ca is not None and cb is not None and ca != cb:
            comm_edges[tuple(sorted((ca, cb)))] += 1

    # community-level layout
    cgraph = ig.Graph.TupleList(
        [(str(a), str(b), float(n)) for (a, b), n in comm_edges.items()],
        directed=False, weights=True,
    )
    clabel = {}
    if cgraph.vcount() > 0:
        # A circular macro-layout keeps the handful of communities readable and
        # never degenerates on a disconnected aggregate graph.
        cl = _rescale(cgraph.layout_circle())
        clabel = {cgraph.vs[i]["name"]: (cl[i][0], cl[i][1])
                  for i in range(cgraph.vcount())}

    comm_path = args.output / "communities.tsv"
    with comm_path.open("w", encoding="utf-8", newline="") as handle:
        header = ["community", "n_families", "n_proteins", "dominant_phylum",
                  "phylum_fraction", "median_length", "cross_phylum_families",
                  "structure_annotated_families", "structure_annotated_fraction",
                  "label", "x", "y"]
        handle.write("\t".join(header) + "\n")
        for row in community_rows:
            x, y = clabel.get(str(row["community"]), ("", ""))
            handle.write(
                f"{row['community']}\t{row['n_families']}\t{row['n_proteins']}\t"
                f"{row['dominant_phylum']}\t{row['phylum_fraction']}\t"
                f"{row['median_length']}\t{row['cross_phylum_families']}\t"
                f"{row['structure_annotated_families']}\t"
                f"{row['structure_annotated_fraction']}\t{row['label']}\t{x}\t{y}\n"
            )

    edges_path = args.output / "community_edges.tsv"
    with edges_path.open("w", encoding="utf-8", newline="") as handle:
        handle.write("community_a\tcommunity_b\tshared_edges\n")
        for (a, b), n in sorted(comm_edges.items(), key=lambda kv: -kv[1]):
            handle.write(f"{a}\t{b}\t{n}\n")

    # Family-level edges restricted to the core set, for the detail plot.
    core_edges_path = args.output / "core_edges.tsv"
    with core_edges_path.open("w", encoding="utf-8", newline="") as handle:
        handle.write("family_a\tfamily_b\tweight\n")
        for a, b, w in family_edges:
            if a in core_set and b in core_set:
                handle.write(f"{a}\t{b}\t{w:.4f}\n")

    summary = {
        "n_families_in_network": len(family_keys),
        "n_edges": len(family_edges),
        "n_communities": len(community_rows),
        "k": args.k,
        "resolution": args.resolution,
        "largest_community": community_rows[0] if community_rows else None,
    }
    (args.output / "network_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"Network written: {len(family_keys)} nodes, {len(family_edges)} edges, "
          f"{len(community_rows)} communities.")
    if community_rows:
        top = community_rows[0]
        print(f"Largest community C{top['community']}: {top['n_proteins']} proteins "
              f"({top['n_families']} families), dominated by {top['dominant_phylum']}.")
    return 0


def _community_label(dominant: str, phylum_frac: float, median_length: float,
                     cross_phylum: int, n_families: int, with_struct: int) -> str:
    bits = []
    if phylum_frac >= 0.6 and dominant:
        bits.append(f"{dominant}-dominated")
    elif dominant:
        bits.append("taxonomically mixed")
    if median_length < 100:
        bits.append("small proteins (<100 aa)")
    elif median_length > 400:
        bits.append("long proteins (>400 aa)")
    if cross_phylum / max(1, n_families) >= 0.5:
        bits.append("broadly cross-phylum conserved")
    if with_struct / max(1, n_families) >= 0.3:
        bits.append("structure-annotated")
    else:
        bits.append("largely unannotated")
    return "; ".join(bits) if bits else "unclassified"


if __name__ == "__main__":
    raise SystemExit(main())
