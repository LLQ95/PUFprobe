#!/usr/bin/env python3
"""Make fair-esm ESMFold importable on Python 3.11 and newer.

fair-esm 2.0.0 declares nested dataclass defaults such as
``structure_module: StructureModuleConfig = StructureModuleConfig()``.
Python 3.11 rejects any non-hashable dataclass default at class-definition
time ("mutable default ... is not allowed: use default_factory"), while the
older check only covered list/dict/set. ESMFold is documented for Python
<=3.9, but these two nested defaults are the sole import-time blockers on
newer Python. Rewriting them as ``field(default_factory=...)`` is
semantically identical and works on every Python version. The patch is
idempotent and safe to re-run.
"""
import os
import sys

PATCHES = {
    os.path.join("esmfold", "v1", "trunk.py"): [
        (
            "from dataclasses import dataclass\n",
            "from dataclasses import dataclass, field\n",
        ),
        (
            "structure_module: StructureModuleConfig = StructureModuleConfig()",
            "structure_module: StructureModuleConfig = field(default_factory=StructureModuleConfig)",
        ),
    ],
    os.path.join("esmfold", "v1", "esmfold.py"): [
        (
            "from dataclasses import dataclass\n",
            "from dataclasses import dataclass, field\n",
        ),
        (
            "trunk: T.Any = FoldingTrunkConfig()",
            "trunk: T.Any = field(default_factory=FoldingTrunkConfig)",
        ),
    ],
}


def patch_file(path, replacements):
    with open(path, "r", encoding="utf-8") as handle:
        text = handle.read()
    original = text
    for old, new in replacements:
        if new in text:
            continue
        if old not in text:
            print("  skip (pattern absent, likely already patched): "
                  + old.splitlines()[0][:70])
            continue
        text = text.replace(old, new)
    if text == original:
        return False
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(text)
    return True


def main():
    try:
        import esm
    except ImportError as error:
        print("ERROR: fair-esm is not installed in this environment: " + str(error))
        return 1
    esm_dir = os.path.dirname(os.path.abspath(esm.__file__))
    changed = []
    for rel, replacements in PATCHES.items():
        target = os.path.join(esm_dir, rel)
        if not os.path.isfile(target):
            print("ERROR: expected fair-esm file not found: " + target)
            return 1
        if patch_file(target, replacements):
            changed.append(target)
            print("patched: " + target)
        else:
            print("already up to date: " + target)
    if changed:
        print("Applied Python 3.11+ dataclass compatibility patch to fair-esm.")
    else:
        print("fair-esm already carries the dataclass compatibility patch.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
