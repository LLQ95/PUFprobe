"""No-op stand-in for ``pytorch_lightning.callbacks.early_stopping``."""

__all__ = ["EarlyStopping"]


class EarlyStopping:
    """No-op base class (training-only). Not instantiated during inference."""
