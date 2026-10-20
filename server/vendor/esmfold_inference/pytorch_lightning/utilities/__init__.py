"""No-op rank-zero logging helpers (training-only, not used by inference)."""

__all__ = ["rank_zero_info", "rank_zero_warn", "rank_zero_only"]


def rank_zero_info(*args, **kwargs):
    pass


def rank_zero_warn(*args, **kwargs):
    pass


def rank_zero_only(fn):
    return fn
