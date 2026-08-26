
from importlib.metadata import version, PackageNotFoundError

__packagename__ = "FastOMA"
try:
    __version__ = version(__packagename__)
except PackageNotFoundError:
    __version__ = "unknown"

import logging
logger = logging.getLogger("FastOMA")
