"""Minimal inference-only stub for NVIDIA ``dllogger``.

Only ``openfold.utils.logger`` (a training callback) imports it. ESMFold
inference never initializes or writes to DLLogger, so the module is inert.
"""

__all__ = [
    "Verbosity",
    "JSONStreamBackend",
    "StdOutBackend",
    "init",
    "log",
    "flush",
]


class Verbosity:
    LOG_ALL = 0
    VERBOSE = 0
    DEFAULT = 1
    INFO = 1
    WARNING = 2
    ERROR = 3


class _BaseBackend:
    def __init__(self, *args, **kwargs):
        pass


class JSONStreamBackend(_BaseBackend):
    pass


class StdOutBackend(_BaseBackend):
    pass


def init(*args, **kwargs):
    pass


def log(*args, **kwargs):
    pass


def flush(*args, **kwargs):
    pass
