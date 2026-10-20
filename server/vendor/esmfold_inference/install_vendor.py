#!/usr/bin/env python3
"""Install the inference-only OpenFold vendor into the active environment.

fair-esm's ESMFold imports OpenFold, but the official OpenFold installer
requires a CUDA toolchain (nvcc) and pulls in training-only packages such as
DeepSpeed and pytorch-lightning. ESMFold inference only uses OpenFold's
pure-PyTorch modules; the single fused CUDA kernel and the training callbacks
are never invoked (inplace_safe=False and the fused attention disabled).

This vendor therefore provides:

  * the OpenFold Python package (Apache-2.0, pinned commit
    4b41059694619831a7db195b7e0988fc4ff3a307), including the
    stereo_chemical_props.txt resource that is downloaded separately by the
    upstream installer;
  * an inert ``attn_core_inplace_cuda`` module that satisfies the
    module-level import without compiling anything;
  * minimal ``pytorch_lightning`` and ``dllogger`` stubs that cover the
    symbols referenced by OpenFold's package initializers.

Installation is decided by the presence of each component in the active
environment's site-packages, not by the import path, so the vendor copy
living next to this script is never mistaken for an installed package. Real
packages are never overwritten, and the operation is idempotent.
"""
import importlib
import os
import shutil
import site
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def site_packages_dir():
    return site.getsitepackages()[0]


def install_pkg(src_dir, dst_dir):
    if os.path.exists(dst_dir):
        shutil.rmtree(dst_dir)
    shutil.copytree(src_dir, dst_dir)


def main():
    site_packages = site_packages_dir()
    actions = []

    # OpenFold: present when its key inference module is installed.
    if not os.path.isfile(
        os.path.join(site_packages, "openfold", "model", "structure_module.py")
    ):
        install_pkg(
            os.path.join(HERE, "openfold"),
            os.path.join(site_packages, "openfold"),
        )
        actions.append("openfold (vendored, pure Python)")

    # Fused CUDA kernel: only stub it when the compiled extension is absent.
    if not os.path.isfile(
        os.path.join(site_packages, "attn_core_inplace_cuda.py")
    ):
        shutil.copy2(
            os.path.join(HERE, "attn_core_inplace_cuda.py"),
            os.path.join(site_packages, "attn_core_inplace_cuda.py"),
        )
        actions.append("attn_core_inplace_cuda (stub)")

    # DLLogger: stub only when the real package is not installed.
    if not os.path.isfile(os.path.join(site_packages, "dllogger.py")):
        shutil.copy2(
            os.path.join(HERE, "dllogger.py"),
            os.path.join(site_packages, "dllogger.py"),
        )
        actions.append("dllogger (stub)")

    # pytorch-lightning: stub only when the real package is not installed.
    if not os.path.isfile(
        os.path.join(site_packages, "pytorch_lightning", "__init__.py")
    ):
        install_pkg(
            os.path.join(HERE, "pytorch_lightning"),
            os.path.join(site_packages, "pytorch_lightning"),
        )
        actions.append("pytorch_lightning (stub)")

    if actions:
        print("Installed into " + site_packages + ":")
        for item in actions:
            print("  - " + item)
    else:
        print("All inference-only OpenFold components already present; nothing to do.")

    # Verify imports resolve from site-packages. Drop this script's directory
    # (sys.path[0]) so the on-disk vendor copy cannot mask a missing install.
    importlib.invalidate_caches()
    saved_path = list(sys.path)
    sys.path[:] = [
        p for p in sys.path
        if p and os.path.abspath(p) != HERE
    ]
    try:
        for name in (
            "openfold.model.structure_module",
            "openfold.utils.kernel.attention_core",
            "attn_core_inplace_cuda",
            "dllogger",
            "pytorch_lightning",
        ):
            importlib.import_module(name)
    finally:
        sys.path[:] = saved_path
    print("Vendor verification: all OpenFold inference imports resolved from site-packages.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
