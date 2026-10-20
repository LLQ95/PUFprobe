"""No-op stand-in for ``pytorch_lightning.utilities.seed`` (training only)."""

__all__ = ["seed_everything"]


def seed_everything(*args, **kwargs):
    return 0
