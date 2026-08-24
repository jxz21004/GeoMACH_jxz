"""Core parametric-geometry classes.

``MACHconfiguration`` adds MPI-based coupling hooks.  Keep it available when
mpi4py is installed, but do not make mpi4py a hard requirement for the base
GeoMACH geometry classes.
"""

try:
    from .MACHconfiguration import MACHconfiguration
except ModuleNotFoundError as err:
    if err.name != "mpi4py":
        raise
    MACHconfiguration = None

from .PGMconfiguration import PGMconfiguration
from .PGMcomponent import PGMcomponent
from .PGMface import PGMface
from .PGMsurf import PGMsurf
from .PGMobject import PGMobject
from .PGMvec import PGMvec
from .PGMproperty import PGMproperty
from .PGMparameter import PGMparameter
from .PGMdv import PGMdv

__all__ = [
    "MACHconfiguration",
    "PGMconfiguration",
    "PGMcomponent",
    "PGMface",
    "PGMsurf",
    "PGMobject",
    "PGMvec",
    "PGMproperty",
    "PGMparameter",
    "PGMdv",
]
