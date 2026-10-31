"""Diagnose GO depth values for the terms that survive the annotate gate."""
import pufprobe
from pufprobe.ontology_graph import load_go_obo

OBO = "/data/liushiwei/PUFprobe/dbs/dbdir/go-basic.obo"
TERMS = [
    "GO:0005886",
    "GO:0009306",
    "GO:0015031",
    "GO:0016020",
    "GO:0019867",
    "GO:0071944",
    "GO:0140352",
]

print("pufprobe file:", pufprobe.__file__)
print("pufprobe version:", pufprobe.__version__)
onto = load_go_obo(OBO)
for term in TERMS:
    primary = onto.resolve(term)
    parents = sorted(onto._parents.get(primary, ()))
    print(
        f"{term}  name={onto.name(term)!r}  depth={onto.depth(term)}  "
        f"parents={parents}"
    )
    for parent in parents:
        print(f"    parent {parent} depth={onto.depth(parent)}")
