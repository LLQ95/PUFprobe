"""Inspect every term and evidence object for one predicted protein."""
import json
from pathlib import Path

WORKDIR = Path("/data/liushiwei/PUFprobe/test_e2e/run")
TARGET = "hp26695@LPJIAC_00556"

pred_path = WORKDIR / "07_rank" / "predictions.json"
preds = json.loads(pred_path.read_text(encoding="utf-8"))

for protein in preds:
    key = f"{protein['sample_id']}@{protein['locus_tag']}"
    if key != TARGET:
        continue
    print("protein", key, "terms:", len(protein["terms"]))
    for term in protein["terms"]:
        print(
            term["term_type"], term["term_id"],
            "name=", repr(term["term_name"]),
            "depth=", term.get("go_depth"),
            "conf=", round(term["confidence"], 3),
        )
        for ev in term["evidence"]:
            details = ev.get("details", {})
            print(
                "   ", ev["channel"],
                "tid=", ev["target_id"],
                "label=", repr(ev["target_label"])[:55],
                "prop=", details.get("propagated", False),
                "xref=", details.get("xref", False),
                "score=", round(ev["score"], 3),
            )
