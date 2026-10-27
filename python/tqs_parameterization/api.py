"""High-level API: Status, result objects, Cache, one-shot functions.

Two tiers, matching the C API. The module-level functions are one-shot: they
compute a shape and keep nothing. :class:`Cache` is the read-only cache on a
primary theta set: build it once, then request any number of shapes against
it. It holds nothing derived from a shape parameter, so no call influences a
later one. Both tiers return the same status and bitwise-identical outputs.

A :class:`Cache` is immutable after construction and may be shared between
threads: ctypes releases the GIL during a call, and the library computes on a
cache concurrently. :meth:`Cache.close` must not race with a compute.

Shape-validation failures come back as result objects carrying a
:class:`Status` with every numeric output zero-filled; :class:`TqsParamError`
is raised only for usage errors — a failed create, a closed handle, or
non-1-D input.

Units and frame: R0 = 1 (volume 4 pi / 3), centre of mass at the origin.

Precondition — every ``params`` and ``thetas`` array must be finite. Release
builds assume no NaN, so non-finite input is undefined behavior.
"""
from __future__ import annotations

import ctypes
from dataclasses import dataclass
from enum import IntEnum
from typing import Optional, Tuple

import numpy as np
import numpy.typing as npt

from ._cdefs import N_PARAMS, N_POINTS_FLOOR, GeometryC, c_dbl_p, configure
from ._libloader import load_library

_lib: Optional[ctypes.CDLL] = None


def _get_lib() -> ctypes.CDLL:
    """Load and configure the shared library once per process."""
    global _lib
    if _lib is None:
        _lib = configure(load_library())
    return _lib


class TqsParamError(RuntimeError):
    """Raised for usage errors: failed create, closed handle, non-1-D input.

    Attributes
    ----------
    status : Status or None
        The library's status code when the error came from a rejected cache
        construction; ``None`` for errors raised on the Python side.
    """

    def __init__(self, message: str, status: "Optional[Status]" = None) -> None:
        super().__init__(message)
        self.status = status


class Status(IntEnum):
    """Status codes of the C API (``TQS_*`` in the header).

    0-5 are the shared shape-parameterization contract codes, identical in
    every library of the family; 100+ are tqs's own.
    """

    valid = 0
    too_many_params = 1
    cache_not_initialized = 2
    invalid_grid = 3
    wrong_param_count = 4
    invalid_init = 5
    invalid_coordinates = 100
    no_geometry = 101
    not_star_convex = 102
    separated = 103


class MidKind(IntEnum):
    """Middle-body kind of a resolved shape (``TQS_MID_*``); ``none`` on failure."""

    none = 0
    spheroid = 1
    cylinder = 2
    one_sheet = 3
    cone = 4
    two_sheet = 5


def status_message(status: int) -> str:
    """Fixed description of a status code.

    Parameters
    ----------
    status : int
        A :class:`Status` value or its integer code.

    Returns
    -------
    str
        The library's static message; ``"unknown status code"`` for codes it
        does not know.
    """
    return _get_lib().tqs_param_status_message(int(status)).decode()


# --- result objects --------------------------------------------------------

class _StatusMixin:
    """``ok`` and ``message`` for every result object."""

    status: Status

    @property
    def ok(self) -> bool:
        """True when the evaluation succeeded."""
        return self.status == Status.valid

    @property
    def message(self) -> str:
        """Description of :attr:`status`."""
        return status_message(self.status)


@dataclass(frozen=True)
class RadiusGridResult(_StatusMixin):
    """R(theta) on a theta set; zero-filled when ``status`` is nonzero."""

    radii: npt.NDArray[np.float64]
    status: Status


@dataclass(frozen=True)
class RadiusDerivativeResult(_StatusMixin):
    """R(theta) and dR/dtheta; both zero-filled on failure."""

    radii: npt.NDArray[np.float64]
    dr_dtheta: npt.NDArray[np.float64]
    status: Status


@dataclass(frozen=True)
class TqsShape(_StatusMixin):
    """Pole radii from the centre of mass, with the R(theta) verdict.

    A pinched or non-star-convex shape is ``not_star_convex`` and a separated
    one ``separated``, with zero poles: :func:`geometry` has the poles and
    flags of every valid geometry.
    """

    r_north: float
    r_south: float
    status: Status


@dataclass(frozen=True)
class RhoZGridResult(_StatusMixin):
    """rho(z) and drho/dz on a uniform grid over [z_south, z_north].

    Both poles are nodes, with rho = drho/dz = 0 there; likewise at a pinch
    point and, in the unchecked form, at every node inside the void of a
    separated shape.
    """

    z: npt.NDArray[np.float64]
    rho: npt.NDArray[np.float64]
    drho_dz: npt.NDArray[np.float64]
    status: Status


@dataclass(frozen=True)
class NeckResult(_StatusMixin):
    """The waist of a middle body with sigma2 < 0 lying strictly inside it.

    ``rho_neck`` is 0 for a pinched or separated shape (the split plane). A
    valid shape without a neck is ``ok`` with ``found`` False and zeros.
    """

    z_neck: float
    rho_neck: float
    found: bool
    status: Status


@dataclass(frozen=True)
class Geometry(_StatusMixin):
    """The resolved geometry.

    Index 0 is the south end spheroid, 1 the north one, 2 the middle body:
    ``rho^2 = a_sq[i] - k[i] (z - l[i])^2``. ``c`` holds the end semi-axes
    along z, ``z_join`` the joins, ``alpha1`` Nix's centre-of-mass coordinate.
    Valid for separated and pinched shapes too.
    """

    u: float
    a_sq: Tuple[float, float, float]
    k: Tuple[float, float, float]
    c: Tuple[float, float]
    l: Tuple[float, float, float]
    z_join: Tuple[float, float]
    r_north: float
    r_south: float
    alpha1: float
    mid_kind: MidKind
    waist_inside: bool
    is_separated: bool
    is_pinched: bool
    status: Status


@dataclass(frozen=True)
class MollerCoordinates(_StatusMixin):
    """Moller's grid coordinates of the shape, as diagnostics.

    ``q2`` is the dimensionless <r^2 P2> over the volume (R0 = 1); Moller's
    charge quadrupole moment is ``Q2 = 2 Z R0^2 q2``. ``d`` is 0 for pinched
    and separated shapes.
    """

    q2: float
    d: float
    alpha_g: float
    eps_f1: float
    eps_f2: float
    status: Status


@dataclass(frozen=True)
class ScissionSigma1(_StatusMixin):
    """The scission map of (sigma2, sigma3, alpha2, alpha3).

    ``sigma1_t`` is Moller's touching-spheroids separation; ``sigma1_star`` is
    where a3^2 = 0 at this sigma2, ``defined`` only for sigma2 < 0 (else 0).
    """

    sigma1_t: float
    sigma1_star: float
    defined: bool
    status: Status


# --- helpers ---------------------------------------------------------------

def theta_grid(n: int) -> npt.NDArray[np.float64]:
    """The open uniform theta grid ``theta_i = i*pi/(n+1)``, ``i = 1..n``.

    Open: neither pole is a node, and no node can round past pi. Pole radii
    come from :func:`shape` / :func:`geometry`.

    Parameters
    ----------
    n : int
        Number of interior points.

    Returns
    -------
    numpy.ndarray
        ``n`` thetas in radians, strictly inside ``(0, pi)``.
    """
    return np.arange(1, int(n) + 1, dtype=np.float64) * np.pi / (int(n) + 1.0)


def _as_1d(values: npt.ArrayLike, name: str) -> npt.NDArray[np.float64]:
    """Contiguous 1-D float64 copy or view of ``values``; TqsParamError otherwise."""
    arr = np.ascontiguousarray(values, dtype=np.float64)
    if arr.ndim != 1:
        raise TqsParamError(f"{name} must be 1-D, got shape {arr.shape}")
    return arr


def _ptr(arr: npt.NDArray[np.float64]) -> ctypes.POINTER(ctypes.c_double):  # type: ignore[valid-type]
    """Raw double pointer to a contiguous float64 array."""
    return arr.ctypes.data_as(c_dbl_p)


def _zeros(n: int) -> npt.NDArray[np.float64]:
    """An output buffer; a negative size gives an empty one and the library judges it."""
    return np.zeros(max(int(n), 0), dtype=np.float64)


def _geometry(g: GeometryC, status: int) -> Geometry:
    return Geometry(u=g.u, a_sq=tuple(g.a_sq), k=tuple(g.k), c=tuple(g.c), l=tuple(g.l),
                    z_join=tuple(g.z_join), r_north=g.r_north, r_south=g.r_south,
                    alpha1=g.alpha1, mid_kind=MidKind(g.mid_kind),
                    waist_inside=bool(g.waist_inside), is_separated=bool(g.is_separated),
                    is_pinched=bool(g.is_pinched), status=Status(status))


# --- one-shot functions ----------------------------------------------------

def radius_grid(params: npt.ArrayLike, thetas: npt.ArrayLike) -> RadiusGridResult:
    """R(theta) for one shape.

    Parameters
    ----------
    params : array_like
        ``[sigma1, sigma2, sigma3, alpha2, alpha3]``, exactly ``N_PARAMS``.
    thetas : array_like
        Evaluation angles in radians, at least one, all inside ``[0, pi]``.

    Returns
    -------
    RadiusGridResult
        Radii from the centre of mass and the status; zero-filled on failure.
    """
    p = _as_1d(params, "params")
    t = _as_1d(thetas, "thetas")
    radii = _zeros(t.size)
    status = ctypes.c_int32(0)
    _get_lib().tqs_param_radius_grid(_ptr(p), p.size, _ptr(t), t.size, _ptr(radii),
                                     ctypes.byref(status))
    return RadiusGridResult(radii=radii, status=Status(status.value))


def radius_and_derivative(params: npt.ArrayLike,
                          thetas: npt.ArrayLike) -> RadiusDerivativeResult:
    """R(theta) and dR/dtheta for one shape.

    Parameters
    ----------
    params : array_like
        ``[sigma1, sigma2, sigma3, alpha2, alpha3]``.
    thetas : array_like
        Evaluation angles in radians, at least one, all inside ``[0, pi]``.

    Returns
    -------
    RadiusDerivativeResult
        Radii, derivatives and status; both zero-filled on failure.
    """
    p = _as_1d(params, "params")
    t = _as_1d(thetas, "thetas")
    radii = _zeros(t.size)
    dr_dtheta = _zeros(t.size)
    status = ctypes.c_int32(0)
    _get_lib().tqs_param_radius_and_derivative(_ptr(p), p.size, _ptr(t), t.size, _ptr(radii),
                                               _ptr(dr_dtheta), ctypes.byref(status))
    return RadiusDerivativeResult(radii=radii, dr_dtheta=dr_dtheta, status=Status(status.value))


def shape(params: npt.ArrayLike) -> TqsShape:
    """Pole radii with the R(theta) verdict (see :class:`TqsShape`)."""
    p = _as_1d(params, "params")
    r_north = ctypes.c_double(0.0)
    r_south = ctypes.c_double(0.0)
    status = ctypes.c_int32(0)
    _get_lib().tqs_param_shape(_ptr(p), p.size, ctypes.byref(r_north), ctypes.byref(r_south),
                               ctypes.byref(status))
    return TqsShape(r_north=r_north.value, r_south=r_south.value, status=Status(status.value))


def _rho_z(params: npt.ArrayLike, n_points: int, checked: bool) -> RhoZGridResult:
    p = _as_1d(params, "params")
    z, rho, drho_dz = _zeros(n_points), _zeros(n_points), _zeros(n_points)
    status = ctypes.c_int32(0)
    entry = _get_lib().tqs_param_rho_z_grid if checked else _get_lib().tqs_param_rho_z_grid_unchecked
    entry(_ptr(p), p.size, int(n_points), _ptr(z), _ptr(rho), _ptr(drho_dz), ctypes.byref(status))
    return RhoZGridResult(z=z, rho=rho, drho_dz=drho_dz, status=Status(status.value))


def rho_z_grid(params: npt.ArrayLike, n_points: int) -> RhoZGridResult:
    """The profile on ``n_points >= N_POINTS_FLOOR`` nodes; ``separated`` for two fragments."""
    return _rho_z(params, n_points, checked=True)


def rho_z_grid_unchecked(params: npt.ArrayLike, n_points: int) -> RhoZGridResult:
    """As :func:`rho_z_grid`, never ``separated``: rho = 0 inside the void."""
    return _rho_z(params, n_points, checked=False)


def neck(params: npt.ArrayLike) -> NeckResult:
    """The analytic neck (see :class:`NeckResult`)."""
    p = _as_1d(params, "params")
    z_neck = ctypes.c_double(0.0)
    rho_neck = ctypes.c_double(0.0)
    found = ctypes.c_int32(0)
    status = ctypes.c_int32(0)
    _get_lib().tqs_param_neck(_ptr(p), p.size, ctypes.byref(z_neck), ctypes.byref(rho_neck),
                              ctypes.byref(found), ctypes.byref(status))
    return NeckResult(z_neck=z_neck.value, rho_neck=rho_neck.value, found=bool(found.value),
                      status=Status(status.value))


def geometry(params: npt.ArrayLike) -> Geometry:
    """The resolved geometry (see :class:`Geometry`)."""
    p = _as_1d(params, "params")
    g = GeometryC()
    status = ctypes.c_int32(0)
    _get_lib().tqs_param_geometry(_ptr(p), p.size, ctypes.byref(g), ctypes.byref(status))
    return _geometry(g, status.value)


def moller_coordinates(params: npt.ArrayLike) -> MollerCoordinates:
    """Moller's grid coordinates (see :class:`MollerCoordinates`)."""
    p = _as_1d(params, "params")
    out = [ctypes.c_double(0.0) for _ in range(5)]
    status = ctypes.c_int32(0)
    _get_lib().tqs_param_moller_coordinates(_ptr(p), p.size, *(ctypes.byref(v) for v in out),
                                            ctypes.byref(status))
    q2, d, alpha_g, eps_f1, eps_f2 = (v.value for v in out)
    return MollerCoordinates(q2=q2, d=d, alpha_g=alpha_g, eps_f1=eps_f1, eps_f2=eps_f2,
                             status=Status(status.value))


def scission_sigma1(params: npt.ArrayLike) -> ScissionSigma1:
    """The scission map (see :class:`ScissionSigma1`); sigma1 is validated but unused."""
    p = _as_1d(params, "params")
    sigma1_t = ctypes.c_double(0.0)
    sigma1_star = ctypes.c_double(0.0)
    defined = ctypes.c_int32(0)
    status = ctypes.c_int32(0)
    _get_lib().tqs_param_scission_sigma1(_ptr(p), p.size, ctypes.byref(sigma1_t),
                                         ctypes.byref(sigma1_star), ctypes.byref(defined),
                                         ctypes.byref(status))
    return ScissionSigma1(sigma1_t=sigma1_t.value, sigma1_star=sigma1_star.value,
                          defined=bool(defined.value), status=Status(status.value))


# --- read-only cache -------------------------------------------------------

class Cache:
    """Read-only cache on a primary theta set.

    Build once, evaluate any number of shapes. Nothing derived from a shape
    parameter is stored, so no call influences a later one. Immutable after
    construction and safe to share between threads; :meth:`close` must not
    race with a compute. The handle is released by :meth:`close`, by leaving
    a ``with`` block, or at garbage collection.

    Parameters
    ----------
    thetas : array_like
        Primary evaluation angles in radians, at least one, all inside
        ``[0, pi]``. See :func:`theta_grid`. The cache keeps its own copy.
    max_params : int
        ``N_PARAMS`` (5): tqs is fixed-arity; the argument exists for
        signature uniformity with the expansion libraries.

    Attributes
    ----------
    max_params, n_thetas : int
        The construction arguments.

    Raises
    ------
    TqsParamError
        If the library refuses the arguments. Its ``status`` carries the code:
        ``too_many_params`` (``max_params > 5``), ``invalid_init``
        (``max_params < 5``) or ``invalid_grid`` (the theta set).
    """

    def __init__(self, thetas: npt.ArrayLike, max_params: int = N_PARAMS) -> None:
        lib = _get_lib()
        t = _as_1d(thetas, "thetas")
        status = ctypes.c_int32(0)
        handle = lib.tqs_param_cache_create(int(max_params), _ptr(t), t.size,
                                            ctypes.byref(status))
        if not handle:
            code = Status(status.value)
            raise TqsParamError(
                f"Cache creation failed (max_params={max_params}, n_thetas={t.size}): "
                f"{status_message(code)} [{code.name}]; need max_params = {N_PARAMS} and "
                f"at least one theta in [0, pi]", status=code)
        self._handle: Optional[ctypes.c_void_p] = ctypes.c_void_p(handle)
        # Bound here so __del__ never needs module globals during shutdown.
        self._destroy = lib.tqs_param_cache_destroy
        self.max_params = int(max_params)
        self.n_thetas = int(t.size)

    def close(self) -> None:
        """Destroy the underlying handle. Idempotent."""
        handle = getattr(self, "_handle", None)
        if handle is not None:
            self._handle = None
            self._destroy(handle)

    def __enter__(self) -> "Cache":
        return self

    def __exit__(self, *exc: object) -> None:
        self.close()

    def __del__(self) -> None:
        self.close()

    def _require_handle(self) -> ctypes.c_void_p:
        handle = getattr(self, "_handle", None)
        if handle is None:
            raise TqsParamError("Cache is closed")
        return handle

    def radius_grid(self, params: npt.ArrayLike) -> RadiusGridResult:
        """R(theta) at the cache's own thetas."""
        handle = self._require_handle()
        p = _as_1d(params, "params")
        radii = _zeros(self.n_thetas)
        status = _get_lib().tqs_param_cache_radius_grid(handle, _ptr(p), p.size, _ptr(radii),
                                                        radii.size)
        return RadiusGridResult(radii=radii, status=Status(status))

    def radius_and_derivative(self, params: npt.ArrayLike) -> RadiusDerivativeResult:
        """R(theta) and dR/dtheta at the cache's own thetas."""
        handle = self._require_handle()
        p = _as_1d(params, "params")
        radii = _zeros(self.n_thetas)
        dr_dtheta = _zeros(self.n_thetas)
        status = _get_lib().tqs_param_cache_radius_and_derivative(
            handle, _ptr(p), p.size, _ptr(radii), _ptr(dr_dtheta), radii.size)
        return RadiusDerivativeResult(radii=radii, dr_dtheta=dr_dtheta, status=Status(status))

    def radius_and_derivative_at_thetas(self, params: npt.ArrayLike,
                                        thetas: npt.ArrayLike) -> RadiusDerivativeResult:
        """R and dR/dtheta at caller thetas (at least one, all inside ``[0, pi]``).

        Given the cache's own thetas this returns the same bits as
        :meth:`radius_and_derivative`.
        """
        handle = self._require_handle()
        p = _as_1d(params, "params")
        t = _as_1d(thetas, "thetas")
        radii = _zeros(t.size)
        dr_dtheta = _zeros(t.size)
        status = _get_lib().tqs_param_cache_radius_and_derivative_at_thetas(
            handle, _ptr(p), p.size, _ptr(t), t.size, _ptr(radii), _ptr(dr_dtheta))
        return RadiusDerivativeResult(radii=radii, dr_dtheta=dr_dtheta, status=Status(status))

    def shape(self, params: npt.ArrayLike) -> TqsShape:
        """Pole radii with the R(theta) verdict."""
        handle = self._require_handle()
        p = _as_1d(params, "params")
        r_north = ctypes.c_double(0.0)
        r_south = ctypes.c_double(0.0)
        status = _get_lib().tqs_param_cache_shape(handle, _ptr(p), p.size,
                                                  ctypes.byref(r_north), ctypes.byref(r_south))
        return TqsShape(r_north=r_north.value, r_south=r_south.value, status=Status(status))

    def _rho_z(self, params: npt.ArrayLike, n_points: int, checked: bool) -> RhoZGridResult:
        handle = self._require_handle()
        p = _as_1d(params, "params")
        z, rho, drho_dz = _zeros(n_points), _zeros(n_points), _zeros(n_points)
        lib = _get_lib()
        entry = lib.tqs_param_cache_rho_z_grid if checked else lib.tqs_param_cache_rho_z_grid_unchecked
        status = entry(handle, _ptr(p), p.size, int(n_points), _ptr(z), _ptr(rho), _ptr(drho_dz))
        return RhoZGridResult(z=z, rho=rho, drho_dz=drho_dz, status=Status(status))

    def rho_z_grid(self, params: npt.ArrayLike, n_points: int) -> RhoZGridResult:
        """The checked profile on ``n_points`` nodes."""
        return self._rho_z(params, n_points, checked=True)

    def rho_z_grid_unchecked(self, params: npt.ArrayLike, n_points: int) -> RhoZGridResult:
        """The unchecked profile on ``n_points`` nodes."""
        return self._rho_z(params, n_points, checked=False)

    def neck(self, params: npt.ArrayLike) -> NeckResult:
        """The analytic neck."""
        handle = self._require_handle()
        p = _as_1d(params, "params")
        z_neck = ctypes.c_double(0.0)
        rho_neck = ctypes.c_double(0.0)
        found = ctypes.c_int32(0)
        status = _get_lib().tqs_param_cache_neck(handle, _ptr(p), p.size, ctypes.byref(z_neck),
                                                 ctypes.byref(rho_neck), ctypes.byref(found))
        return NeckResult(z_neck=z_neck.value, rho_neck=rho_neck.value, found=bool(found.value),
                          status=Status(status))

    def geometry(self, params: npt.ArrayLike) -> Geometry:
        """The resolved geometry."""
        handle = self._require_handle()
        p = _as_1d(params, "params")
        g = GeometryC()
        status = _get_lib().tqs_param_cache_geometry(handle, _ptr(p), p.size, ctypes.byref(g))
        return _geometry(g, status)

    def moller_coordinates(self, params: npt.ArrayLike) -> MollerCoordinates:
        """Moller's grid coordinates."""
        handle = self._require_handle()
        p = _as_1d(params, "params")
        out = [ctypes.c_double(0.0) for _ in range(5)]
        status = _get_lib().tqs_param_cache_moller_coordinates(
            handle, _ptr(p), p.size, *(ctypes.byref(v) for v in out))
        q2, d, alpha_g, eps_f1, eps_f2 = (v.value for v in out)
        return MollerCoordinates(q2=q2, d=d, alpha_g=alpha_g, eps_f1=eps_f1, eps_f2=eps_f2,
                                 status=Status(status))

    def scission_sigma1(self, params: npt.ArrayLike) -> ScissionSigma1:
        """The scission map."""
        handle = self._require_handle()
        p = _as_1d(params, "params")
        sigma1_t = ctypes.c_double(0.0)
        sigma1_star = ctypes.c_double(0.0)
        defined = ctypes.c_int32(0)
        status = _get_lib().tqs_param_cache_scission_sigma1(
            handle, _ptr(p), p.size, ctypes.byref(sigma1_t), ctypes.byref(sigma1_star),
            ctypes.byref(defined))
        return ScissionSigma1(sigma1_t=sigma1_t.value, sigma1_star=sigma1_star.value,
                              defined=bool(defined.value), status=Status(status))


__all__ = [
    "Cache", "Geometry", "MidKind", "MollerCoordinates", "NeckResult",
    "RadiusDerivativeResult", "RadiusGridResult", "RhoZGridResult", "ScissionSigma1",
    "Status", "TqsParamError", "TqsShape", "N_PARAMS", "N_POINTS_FLOOR",
    "geometry", "moller_coordinates", "neck", "radius_and_derivative", "radius_grid",
    "rho_z_grid", "rho_z_grid_unchecked", "scission_sigma1", "shape", "status_message",
    "theta_grid",
]
