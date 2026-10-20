# Inference-only OpenFold vendor for ESMFold

This directory carries the small set of third-party modules that fair-esm's
ESMFold needs to load and run, prepared so that **no CUDA compiler and no
training-only dependencies are required**.

## Why this is vendored

`pip install fair-esm[esmfold]` pulls OpenFold from a specific Git commit and
builds a fused CUDA extension (`attn_core_inplace_cuda`) with `nvcc`, and it
also pulls DeepSpeed, pytorch-lightning and other training packages. On a
deployment server that only runs ESMFold inference, this is both unnecessary
and fragile (the toolchain is often absent, and the pinned DeepSpeed
0.5.9 conflicts with modern PyTorch).

During ESMFold inference the fused kernel and the training callbacks are never
executed: OpenFold runs with `inplace_safe=False` and
`use_memory_efficient_kernel=False`, taking the pure-PyTorch attention and
softmax branches. The only hard requirement is that the corresponding modules
can be imported.

## Contents

| Path | Origin | Purpose |
|---|---|---|
| `openfold/` | OpenFold commit `4b41059694619831a7db195b7e0988fc4ff3a307` (Apache-2.0) | Pure-PyTorch OpenFold modules, including `resources/stereo_chemical_props.txt` (sourced from OpenStructure, as the upstream installer does) and `LICENSE`. |
| `attn_core_inplace_cuda.py` | PUFprobe stub | Satisfies the module-level import; raises if the training-only kernel is ever invoked. |
| `pytorch_lightning/` | PUFprobe stub | The handful of symbols (`Callback`, `EarlyStopping`, `rank_zero_info`, `seed_everything`) referenced by OpenFold initializers. |
| `dllogger.py` | PUFprobe stub | Inert replacement for NVIDIA DLLogger. |
| `install_vendor.py` | PUFprobe | Copies the components into the active `site-packages`, skipping any that are already installed. |

## Install

Run from the PUFprobe root, with the target environment active:

```bash
python server/vendor/esmfold_inference/install_vendor.py
```

This is normally invoked by `server/install_esmfold_inference.sh`, which also
installs the remaining pure-Python dependencies and applies the fair-esm
Python 3.11+ dataclass patch. To train OpenFold, install the real OpenFold
toolchain instead; the installer never overwrites packages that already import
successfully.

## License

OpenFold is distributed under the Apache License 2.0; its `LICENSE` is kept
inside `openfold/`. The three stub modules are part of PUFprobe (MIT).
