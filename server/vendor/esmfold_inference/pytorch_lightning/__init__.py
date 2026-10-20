"""Minimal inference-only stub for ``pytorch_lightning``.

OpenFold's package initializers import its training-only helper modules
(``openfold/utils/callbacks.py``, ``utils/logger.py`` and ``utils/seed.py``),
which reference pytorch_lightning. ESMFold inference never instantiates these
callbacks or seeding helpers, so this stub exposes only the surface that makes
those modules import. Install the real pytorch_lightning package if you need
to train OpenFold.
"""

__all__ = ["Callback", "EarlyStopping"]


class Callback:
    """No-op stand-in for ``pytorch_lightning.Callback`` (training only)."""


from pytorch_lightning.callbacks.early_stopping import EarlyStopping  # noqa: E402
