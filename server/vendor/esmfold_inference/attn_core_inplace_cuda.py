"""Pure-Python stub for OpenFold's compiled ``attn_core_inplace_cuda`` extension.

OpenFold normally builds this fused in-place attention kernel from
``openfold/utils/kernel/csrc/`` with nvcc at install time. The kernel is only
used on training paths: the fused memory-efficient attention when
``use_memory_efficient_kernel=True`` and the in-place softmax when
``inplace_safe=True``. ESMFold runs inference with both disabled, so OpenFold
takes the pure-PyTorch ``_attention`` / ``softmax`` branches and never invokes
these functions.

This module exists solely so that the unconditional module-level imports in
``openfold.model.structure_module`` and
``openfold.utils.kernel.attention_core`` succeed without a CUDA toolchain.
"""

__all__ = ["forward_", "backward_"]


def _raise(*args, **kwargs):
    raise RuntimeError(
        "attn_core_inplace_cuda is a training-only CUDA kernel that is not "
        "compiled in this inference-only OpenFold installation. It should "
        "never be called during ESMFold inference (inplace_safe=False and "
        "use_memory_efficient_kernel=False). If you reach this error, pass "
        "the pure-PyTorch options instead of the fused kernel."
    )


def forward_(*args, **kwargs):
    _raise()


def backward_(*args, **kwargs):
    _raise()
