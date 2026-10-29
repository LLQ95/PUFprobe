#!/usr/bin/env python
"""Audit which evidence channels are present in one or more search files.

A learning-to-rank model can only use channels that were present when it was
trained. When the training set and the runtime data contain different
channels (for example plm_neighbor is present at runtime but absent in
training), the model silently ignores the extra channel. This script reports
the channel and evidence-group composition of each file so that mismatch is
visible before training or deployment.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

# Evidence channel -> independent evidence group (mirrors models.CHANNEL_GROUPS)
CHANNEL_GROUPS = {
    "plm_neighbor": "sequence",
    "foldseek": "structure",
    "remote_homology": "sequence_profile",
    "pocket": "structure",
    "context": "context",
    "database": "database",
    "eggnog": "database",
}


def audit(path: Path) -> None:
    with path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)

    channels: Counter = Counter()
    groups: Counter = Counter()
    proteins_with = Counter()
    for evs in data.values():
        seen_ch = set()
        for ev in evs:
            ch = ev.get("channel", "unknown")
            channels[ch] += 1
            groups[CHANNEL_GROUPS.get(ch, "unknown")] += 1
            seen_ch.add(ch)
        for ch in seen_ch:
            proteins_with[ch] += 1

    print(f"=== {path} ===")
    print(f"proteins: {len(data)}")
    print("evidence rows by channel:")
    for ch, n in sorted(channels.items()):
        print(f"  {ch:16s} {n:>9d} rows  ({proteins_with[ch]} proteins)")
    print("evidence rows by independent group:")
    for grp, n in sorted(groups.items()):
        print(f"  {grp:16s} {n:>9d} rows")
    print()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("evidence", nargs="+", type=Path)
    args = parser.parse_args()
    for path in args.evidence:
        audit(path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
