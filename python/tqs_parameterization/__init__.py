"""Python bindings for the three-quadratic-surface (3QS) nuclear-shape parameterization.

Nix's coordinates ``params = [sigma1, sigma2, sigma3, alpha2, alpha3]``; the
sphere is ``[0, 1, 1, 0, 0]``. Two tiers: one-shot module-level functions,
and :class:`Cache`, a read-only cache on a primary theta set, built once and
shared freely between threads. Validation failures come back as a
:class:`Status` on the result object; only usage errors raise
:class:`TqsParamError`.
"""
from ._cdefs import N_PARAMS, N_POINTS_FLOOR
from ._libloader import load_library
from .api import (
    Cache,
    Geometry,
    MidKind,
    MollerCoordinates,
    NeckResult,
    RadiusDerivativeResult,
    RadiusGridResult,
    RhoZGridResult,
    ScissionSigma1,
    Status,
    TqsParamError,
    TqsShape,
    geometry,
    moller_coordinates,
    neck,
    radius_and_derivative,
    radius_grid,
    rho_z_grid,
    rho_z_grid_unchecked,
    scission_sigma1,
    shape,
    status_message,
    theta_grid,
)

__all__ = [
    "Cache", "Geometry", "MidKind", "MollerCoordinates", "NeckResult",
    "RadiusDerivativeResult", "RadiusGridResult", "RhoZGridResult", "ScissionSigma1",
    "Status", "TqsParamError", "TqsShape",
    "geometry", "moller_coordinates", "neck", "radius_and_derivative", "radius_grid",
    "rho_z_grid", "rho_z_grid_unchecked", "scission_sigma1", "shape",
    "status_message", "theta_grid",
    "N_PARAMS", "N_POINTS_FLOOR",
    "load_library",
]
