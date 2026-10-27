"""ctypes signatures for the tqs-parameterization 0.1.0 C API.

One entry per exported symbol in ``include/tqs_parameterization.h``. Every
integer the C API takes or returns is ``int32_t``; every buffer is ``double*``.
"""
from __future__ import annotations

import ctypes

c_dbl_p = ctypes.POINTER(ctypes.c_double)
c_i32 = ctypes.c_int32
c_i32_p = ctypes.POINTER(ctypes.c_int32)
c_void = ctypes.c_void_p

#: The fixed parameter count, ``[sigma1, sigma2, sigma3, alpha2, alpha3]``
#: (``TQS_PARAM_N_PARAMS``).
N_PARAMS = 5
#: Fewest nodes of a rho(z) grid (``TQS_PARAM_N_POINTS_FLOOR``).
N_POINTS_FLOOR = 3


class GeometryC(ctypes.Structure):
    """``tqs_param_geometry_t``, field for field in header order."""

    _fields_ = [
        ("u", ctypes.c_double),
        ("a_sq", ctypes.c_double * 3),
        ("k", ctypes.c_double * 3),
        ("c", ctypes.c_double * 2),
        ("l", ctypes.c_double * 3),
        ("z_join", ctypes.c_double * 2),
        ("r_north", ctypes.c_double),
        ("r_south", ctypes.c_double),
        ("alpha1", ctypes.c_double),
        ("mid_kind", ctypes.c_int32),
        ("waist_inside", ctypes.c_int32),
        ("is_separated", ctypes.c_int32),
        ("is_pinched", ctypes.c_int32),
    ]


c_geom_p = ctypes.POINTER(GeometryC)

# name -> (argtypes, restype)
_SIGNATURES = {
    # --- diagnostics ---
    "tqs_param_status_message": ([c_i32], ctypes.c_char_p),

    # --- cache lifecycle; the trailing int32_t* receives the create status ---
    "tqs_param_cache_create": ([c_i32, c_dbl_p, c_i32, c_i32_p], c_void),
    "tqs_param_cache_destroy": ([c_void], None),

    # --- cached computes, status is the return value ---
    "tqs_param_cache_radius_grid": ([c_void, c_dbl_p, c_i32, c_dbl_p, c_i32], c_i32),
    "tqs_param_cache_radius_and_derivative": (
        [c_void, c_dbl_p, c_i32, c_dbl_p, c_dbl_p, c_i32], c_i32),
    "tqs_param_cache_radius_and_derivative_at_thetas": (
        [c_void, c_dbl_p, c_i32, c_dbl_p, c_i32, c_dbl_p, c_dbl_p], c_i32),
    "tqs_param_cache_shape": ([c_void, c_dbl_p, c_i32, c_dbl_p, c_dbl_p], c_i32),
    "tqs_param_cache_rho_z_grid": (
        [c_void, c_dbl_p, c_i32, c_i32, c_dbl_p, c_dbl_p, c_dbl_p], c_i32),
    "tqs_param_cache_rho_z_grid_unchecked": (
        [c_void, c_dbl_p, c_i32, c_i32, c_dbl_p, c_dbl_p, c_dbl_p], c_i32),
    "tqs_param_cache_neck": ([c_void, c_dbl_p, c_i32, c_dbl_p, c_dbl_p, c_i32_p], c_i32),
    "tqs_param_cache_geometry": ([c_void, c_dbl_p, c_i32, c_geom_p], c_i32),
    "tqs_param_cache_moller_coordinates": (
        [c_void, c_dbl_p, c_i32, c_dbl_p, c_dbl_p, c_dbl_p, c_dbl_p, c_dbl_p], c_i32),
    "tqs_param_cache_scission_sigma1": (
        [c_void, c_dbl_p, c_i32, c_dbl_p, c_dbl_p, c_i32_p], c_i32),

    # --- one-shot computes, trailing nullable status out-parameter ---
    "tqs_param_radius_grid": ([c_dbl_p, c_i32, c_dbl_p, c_i32, c_dbl_p, c_i32_p], None),
    "tqs_param_radius_and_derivative": (
        [c_dbl_p, c_i32, c_dbl_p, c_i32, c_dbl_p, c_dbl_p, c_i32_p], None),
    "tqs_param_shape": ([c_dbl_p, c_i32, c_dbl_p, c_dbl_p, c_i32_p], None),
    "tqs_param_rho_z_grid": ([c_dbl_p, c_i32, c_i32, c_dbl_p, c_dbl_p, c_dbl_p, c_i32_p], None),
    "tqs_param_rho_z_grid_unchecked": (
        [c_dbl_p, c_i32, c_i32, c_dbl_p, c_dbl_p, c_dbl_p, c_i32_p], None),
    "tqs_param_neck": ([c_dbl_p, c_i32, c_dbl_p, c_dbl_p, c_i32_p, c_i32_p], None),
    "tqs_param_geometry": ([c_dbl_p, c_i32, c_geom_p, c_i32_p], None),
    "tqs_param_moller_coordinates": (
        [c_dbl_p, c_i32, c_dbl_p, c_dbl_p, c_dbl_p, c_dbl_p, c_dbl_p, c_i32_p], None),
    "tqs_param_scission_sigma1": ([c_dbl_p, c_i32, c_dbl_p, c_dbl_p, c_i32_p, c_i32_p], None),
}


def configure(lib: ctypes.CDLL) -> ctypes.CDLL:
    """Set argtypes/restype on every exported symbol (idempotent).

    Parameters
    ----------
    lib : ctypes.CDLL
        The loaded shared library.

    Returns
    -------
    ctypes.CDLL
        The same object, with every signature attached.

    Raises
    ------
    AttributeError
        If the library lacks a 0.1.0 symbol.
    """
    for name, (argtypes, restype) in _SIGNATURES.items():
        func = getattr(lib, name)
        func.argtypes = argtypes
        func.restype = restype
    return lib
