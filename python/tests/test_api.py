"""Contract tests for the 0.1.0 Python surface.

Shape-validation failures are results, not exceptions: the buffers come back
zero-filled with a nonzero :class:`Status`. Only usage errors — a failed
create, a closed handle, non-1-D input — raise ``TqsParamError``.

The golden block is captured by ``tests/golden_capture.f08`` through the same
C entries these tests call. The recapture procedure lives in that file's
header; never loosen ``GOLDEN_ATOL`` to make a recapture fit.
"""
from __future__ import annotations

import math
import threading

import numpy as np
import pytest

import tqs_parameterization as tp

# --- fixtures (params = [sigma1, sigma2, sigma3, alpha2, alpha3]) ----------

SPHERE = [0.0, 1.0, 1.0, 0.0, 0.0]
SPHERE_B = [0.0, -1.0, 1.0, 0.0, 0.0]
SYM_NECKED = [1.6, -0.6, 1.0, 0.0, 0.0]
ASYM_NECKED = [1.6, -0.6, 1.0, 0.3, 0.2]
SPHEROID_MID = [0.8, 0.5, 1.0, 0.2, 0.1]
PEAR = [1.0, -1.0, 1.0, -1.4, 0.0]
PINCHED = [2.0 * math.sqrt(2.0), -1.0, 1.0, 0.0, 0.0]
SEPARATED = [3.0, -1.0, 1.0, 0.0, 0.0]
STAR_FAIL = [3.0, -0.6, 1.0, 0.9, 0.0]
BAD_K = [0.0, 1.0, 0.5, 0.0, 1.0]
NO_GEOM = [0.3, 1.5, 1.5, 0.0, 0.0]

N_Z = 181
GOLDEN_IDX = [0, 45, 90, 135, 180]
# pi/8, pi/2, 7pi/8 as pre-rounded literals, bit-identical to golden_capture.f08.
DERIV_THETAS = [0.39269908169872414, 1.5707963267948966, 2.748893571891069]
# Closed form throughout: the Debug-capture vs shipped-build gap is a few ulp
# of values of order unity, two orders of magnitude inside this.
GOLDEN_ATOL = 1e-13


def _same_bits(a: np.ndarray, b: np.ndarray) -> bool:
    """True when ``a`` and ``b`` hold the same IEEE754 bit patterns.

    ``np.array_equal`` calls ``+0.0`` and ``-0.0`` equal; the contract's
    bitwise claims need the bytes.
    """
    a, b = np.asarray(a), np.asarray(b)
    return a.shape == b.shape and a.dtype == b.dtype and a.tobytes() == b.tobytes()


def _inclusive_grid(n: int) -> np.ndarray:
    """``n`` uniform nodes on the CLOSED ``[0, pi]``, both endpoints pinned."""
    grid = np.arange(n, dtype=np.float64) * np.pi / (n - 1.0)
    grid[0] = 0.0
    grid[-1] = np.pi
    return grid


# --- golden block (tests/golden_capture.f08, Debug) ------------------------
SYM_NECKED_MID_KIND = 3
SYM_NECKED_GEOMETRY = [
     7.8796398791780453E-001,
    -3.9398199395890227E-001,
     3.9398199395890227E-001,
     0.0000000000000000E+000,
     4.7187430715405076E-001,
     1.4183351782520481E+000,
     1.4183351782520481E+000,
     0.0000000000000000E+000,
]
SYM_NECKED_NECK_FOUND = True
SYM_NECKED_NECK = [
     0.0000000000000000E+000,
     6.8693107889660276E-001,
]
SYM_NECKED_MOLLER = [
     3.8054859747368652E-001,
     1.3738621577932055E+000,
     0.0000000000000000E+000,
     0.0000000000000000E+000,
     0.0000000000000000E+000,
]
SYM_NECKED_SCISSION = [
     2.0000000000000000E+000,
     3.2659863237109041E+000,
]
SYM_NECKED_DERIV_R = [
     1.3325164119773067E+000,
     6.8693107889660288E-001,
     1.3325164119773067E+000,
]
SYM_NECKED_DERIV_DR = [
    -4.2852129945887241E-001,
    -6.7299835760445010E-017,
     4.2852129945887252E-001,
]
SYM_NECKED_Z = [
    -1.4183351782520481E+000,
    -7.0916758912602407E-001,
     0.0000000000000000E+000,
     7.0916758912602385E-001,
     1.4183351782520481E+000,
]
SYM_NECKED_RHO = [
     0.0000000000000000E+000,
     7.8401426887064796E-001,
     6.8693107889660276E-001,
     7.8401426887064807E-001,
     0.0000000000000000E+000,
]
ASYM_NECKED_MID_KIND = 3
ASYM_NECKED_GEOMETRY = [
     7.8710105860136714E-001,
    -2.3669474449950981E-001,
     5.5722915207197676E-001,
     2.9002485982038290E-001,
     4.5520073118600868E-001,
     1.5002900817541032E+000,
     1.3287872784416646E+000,
     1.3427056014010980E-001,
]
ASYM_NECKED_NECK_FOUND = True
ASYM_NECKED_NECK = [
     2.9002485982038290E-001,
     6.7468565360915200E-001,
]
ASYM_NECKED_MOLLER = [
     3.6491970920562561E-001,
     1.3493713072183040E+000,
     1.7473143125602089E-001,
    -4.8027459838401433E-002,
     5.2209779853347064E-002,
]
ASYM_NECKED_SCISSION = [
     1.9943000320986402E+000,
     3.2576721506375197E+000,
]
ASYM_NECKED_DERIV_R = [
     1.3689000198337100E+000,
     7.1110434028650704E-001,
     1.2773307541680223E+000,
]
ASYM_NECKED_DERIV_DR = [
    -6.4462999699928147E-001,
     1.7401491589222964E-001,
     2.6242140503223604E-001,
]
ASYM_NECKED_Z = [
    -1.3287872784416646E+000,
    -6.2151793839272262E-001,
     8.5751401656219306E-002,
     7.9302074170516113E-001,
     1.5002900817541032E+000,
]
ASYM_NECKED_RHO = [
     0.0000000000000000E+000,
     8.3785188489471152E-001,
     6.9299157181903481E-001,
     7.2360703070742760E-001,
     0.0000000000000000E+000,
]
PEAR_MID_KIND = 5
PEAR_GEOMETRY = [
     7.6189785796221887E-001,
    -1.0852386959656890E+000,
    -7.0428976698457957E-001,
    -1.4280927320486874E+000,
    -6.0951276326578763E-002,
     1.0129064677087372E+000,
     1.1596933032034649E+000,
    -4.7438869544571427E-001,
]
PEAR_NECK_FOUND = False
PEAR_NECK = [
     0.0000000000000000E+000,
     0.0000000000000000E+000,
]
PEAR_MOLLER = [
     1.8018701068839349E-002,
     4.7580505979580862E-001,
    -8.6196772416877399E-001,
     0.0000000000000000E+000,
     0.0000000000000000E+000,
]
PEAR_SCISSION = [
     1.8515630385456960E+000,
     2.6185055607000609E+000,
]
PEAR_DERIV_R = [
     1.0113930472615786E+000,
     9.9320160251850120E-001,
     1.0391289842351750E+000,
]
PEAR_DERIV_DR = [
    -7.6028974854783857E-003,
    -1.9513198079528338E-002,
     3.3464067012379567E-001,
]
PEAR_Z = [
    -1.1596933032034649E+000,
    -6.1654336047541436E-001,
    -7.3393417747363854E-002,
     4.6975652498068676E-001,
     1.0129064677087372E+000,
]
PEAR_RHO = [
     0.0000000000000000E+000,
     7.6306110008560424E-001,
     9.8903920492577246E-001,
     8.8550049956258425E-001,
     0.0000000000000000E+000,
]
SPHEROID_MID_MID_KIND = 1
SPHEROID_MID_GEOMETRY = [
     8.3601775350464269E-001,
    -5.0348624297648348E-001,
     8.6006922601378832E-001,
    -4.6242079000453747E-002,
     8.2357526403033798E-001,
     1.1967832043907167E+000,
     1.1414423779711795E+000,
     5.8199929846543777E-002,
]
SPHEROID_MID_NECK_FOUND = False
SPHEROID_MID_NECK = [
     0.0000000000000000E+000,
     0.0000000000000000E+000,
]
SPHEROID_MID_MOLLER = [
     1.1780433471663496E-001,
     1.2851072690355909E+000,
     1.2482767431630462E-001,
    -2.4493454021572478E-002,
     2.5536094974404153E-002,
]
SPHEROID_MID_SCISSION = [
     1.9968611582228932E+000,
     0.0000000000000000E+000,
]
SPHEROID_MID_DERIV_R = [
     1.1479408018068484E+000,
     9.0692121989464769E-001,
     1.1177928905954149E+000,
]
SPHEROID_MID_DERIV_DR = [
    -2.4013420570677610E-001,
     2.3121039500226846E-002,
     1.1961876746044343E-001,
]
SPHEROID_MID_Z = [
    -1.1414423779711795E+000,
    -5.5688598238070552E-001,
     2.7670413209768485E-002,
     6.1222680880024249E-001,
     1.1967832043907167E+000,
]
SPHEROID_MID_RHO = [
     0.0000000000000000E+000,
     8.3164194417308379E-001,
     9.0600426918308441E-001,
     7.7896381554575111E-001,
     0.0000000000000000E+000,
]

# --- end golden block ------------------------------------------------------

GOLDEN_SHAPES = {
    "SYM_NECKED": SYM_NECKED,
    "ASYM_NECKED": ASYM_NECKED,
    "PEAR": PEAR,
    "SPHEROID_MID": SPHEROID_MID,
}


# --- theta grid ------------------------------------------------------------

def test_theta_grid_is_open_and_uniform() -> None:
    thetas = tp.theta_grid(64)
    assert thetas.size == 64
    assert thetas[0] > 0.0 and thetas[-1] < np.pi
    assert np.allclose(np.diff(thetas), np.pi / 65.0, rtol=0.0, atol=1e-14)


@pytest.mark.parametrize("n", [1, 2, 7, 181, 2001, 100001])
def test_theta_grid_never_leaves_the_domain(n: int) -> None:
    thetas = tp.theta_grid(n)
    assert thetas.size == n
    assert np.all(thetas > 0.0) and np.all(thetas < np.pi)


def test_theta_grid_is_accepted_by_cache() -> None:
    with tp.Cache(tp.theta_grid(2001)) as cache:
        assert cache.n_thetas == 2001
        assert cache.radius_grid(SYM_NECKED).ok


# --- happy path ------------------------------------------------------------

def test_sphere_is_unit() -> None:
    grid = tp.theta_grid(64)
    res = tp.radius_grid(SPHERE, grid)
    assert res.ok
    assert np.allclose(res.radii, 1.0, rtol=0.0, atol=1e-14)
    shp = tp.shape(SPHERE)
    assert shp.ok and abs(shp.r_north - 1.0) <= 1e-14 and abs(shp.r_south - 1.0) <= 1e-14
    nk = tp.neck(SPHERE)
    assert nk.ok and not nk.found
    geo = tp.geometry(SPHERE)
    assert geo.ok and abs(geo.u - 1.0) <= 1e-14
    res_b = tp.radius_grid(SPHERE_B, grid)
    assert res_b.ok
    assert np.allclose(res_b.radii, res.radii, rtol=0.0, atol=1e-15)


@pytest.mark.parametrize("name", list(GOLDEN_SHAPES))
def test_goldens_match_fortran(name: str) -> None:
    params = GOLDEN_SHAPES[name]
    g = globals()
    geo = tp.geometry(params)
    assert geo.ok
    assert int(geo.mid_kind) == g[f"{name}_MID_KIND"]
    got = [geo.u, geo.z_join[0], geo.z_join[1], geo.l[2], geo.a_sq[2], geo.r_north,
           geo.r_south, geo.alpha1]
    assert np.allclose(got, g[f"{name}_GEOMETRY"], rtol=0.0, atol=GOLDEN_ATOL)

    nk = tp.neck(params)
    assert nk.ok and nk.found == g[f"{name}_NECK_FOUND"]
    assert np.allclose([nk.z_neck, nk.rho_neck], g[f"{name}_NECK"], rtol=0.0, atol=GOLDEN_ATOL)

    mo = tp.moller_coordinates(params)
    assert mo.ok
    assert np.allclose([mo.q2, mo.d, mo.alpha_g, mo.eps_f1, mo.eps_f2], g[f"{name}_MOLLER"],
                       rtol=0.0, atol=GOLDEN_ATOL)

    sc = tp.scission_sigma1(params)
    assert sc.ok
    assert np.allclose([sc.sigma1_t, sc.sigma1_star], g[f"{name}_SCISSION"], rtol=0.0, atol=GOLDEN_ATOL)

    rd = tp.radius_and_derivative(params, DERIV_THETAS)
    assert rd.ok
    assert np.allclose(rd.radii, g[f"{name}_DERIV_R"], rtol=0.0, atol=GOLDEN_ATOL)
    assert np.allclose(rd.dr_dtheta, g[f"{name}_DERIV_DR"], rtol=0.0, atol=GOLDEN_ATOL)

    prof = tp.rho_z_grid(params, N_Z)
    assert prof.ok
    assert np.allclose(prof.z[GOLDEN_IDX], g[f"{name}_Z"], rtol=0.0, atol=GOLDEN_ATOL)
    assert np.allclose(prof.rho[GOLDEN_IDX], g[f"{name}_RHO"], rtol=0.0, atol=GOLDEN_ATOL)


def test_derivative_matches_finite_difference() -> None:
    thetas = tp.theta_grid(90)
    h = 1e-6
    rd = tp.radius_and_derivative(ASYM_NECKED, thetas)
    plus = tp.radius_grid(ASYM_NECKED, thetas + h)
    minus = tp.radius_grid(ASYM_NECKED, thetas - h)
    assert rd.ok and plus.ok and minus.ok
    fd = (plus.radii - minus.radii) / (2.0 * h)
    assert np.all(np.abs(rd.dr_dtheta - fd) <= 1e-6 * (1.0 + np.abs(rd.dr_dtheta)))


@pytest.mark.parametrize("params", [SYM_NECKED, ASYM_NECKED, PEAR, SPHEROID_MID])
def test_profile_tips_are_exact_zero(params) -> None:
    prof = tp.rho_z_grid(params, 101)
    assert prof.ok
    assert prof.rho[0] == 0.0 and prof.rho[-1] == 0.0
    assert prof.drho_dz[0] == 0.0 and prof.drho_dz[-1] == 0.0
    assert np.all(prof.rho[1:-1] > 0.0)


def test_unchecked_profile_draws_a_separated_shape() -> None:
    checked = tp.rho_z_grid(SEPARATED, 101)
    assert checked.status == tp.Status.separated
    raw = tp.rho_z_grid_unchecked(SEPARATED, 101)
    assert raw.ok
    assert np.any(raw.rho[1:-1] == 0.0)
    assert raw.rho[1] > 0.0 and raw.rho[-2] > 0.0


def test_unchecked_equals_checked_on_connected() -> None:
    a = tp.rho_z_grid(ASYM_NECKED, 257)
    b = tp.rho_z_grid_unchecked(ASYM_NECKED, 257)
    assert a.ok and b.ok
    assert _same_bits(a.z, b.z) and _same_bits(a.rho, b.rho) and _same_bits(a.drho_dz, b.drho_dz)


def test_geometry_flags() -> None:
    sep = tp.geometry(SEPARATED)
    assert sep.ok and sep.is_separated and sep.r_north > 0.0 and sep.r_south > 0.0
    pin = tp.geometry(PINCHED)
    assert pin.ok and pin.is_pinched and pin.mid_kind == tp.MidKind.cone
    pear = tp.geometry(PEAR)
    assert pear.ok and pear.mid_kind == tp.MidKind.two_sheet and not pear.waist_inside
    sym = tp.geometry(SYM_NECKED)
    assert sym.ok and sym.mid_kind == tp.MidKind.one_sheet and sym.waist_inside
    assert isinstance(sym.a_sq, tuple) and len(sym.a_sq) == 3 and len(sym.c) == 2
    bad = tp.geometry(NO_GEOM)
    assert bad.status == tp.Status.no_geometry and bad.mid_kind == tp.MidKind.none


def test_scission_sigma1() -> None:
    sc = tp.scission_sigma1([1.0, -1.0, 1.0, 0.0, 0.0])
    assert sc.ok and sc.defined
    assert abs(sc.sigma1_t - 2.0) <= 1e-14
    assert abs(sc.sigma1_star - 2.0 * math.sqrt(2.0)) <= 1e-14
    for s2 in (0.0, 0.5):
        sc = tp.scission_sigma1([1.0, s2, 1.0, 0.0, 0.0])
        assert sc.ok and not sc.defined and sc.sigma1_star == 0.0 and sc.sigma1_t > 0.0


def test_star_fail_is_a_status() -> None:
    assert tp.radius_grid(STAR_FAIL, tp.theta_grid(16)).status == tp.Status.not_star_convex
    assert tp.geometry(STAR_FAIL).ok
    nk = tp.neck(STAR_FAIL)
    assert nk.ok and nk.found and nk.rho_neck > 0.0


def _all_one_shot(params, thetas=None, n_points: int = 33):
    """Every one-shot output for ``params``, in a fixed order."""
    thetas = tp.theta_grid(8) if thetas is None else thetas
    return [
        tp.radius_grid(params, thetas),
        tp.radius_and_derivative(params, thetas),
        tp.shape(params),
        tp.rho_z_grid(params, n_points),
        tp.rho_z_grid_unchecked(params, n_points),
        tp.neck(params),
        tp.geometry(params),
        tp.moller_coordinates(params),
        tp.scission_sigma1(params),
    ]


def _all_cached(cache, params, n_points: int = 33):
    """Every cached output for ``params``, in the one-shot order plus at-thetas."""
    return [
        cache.radius_grid(params),
        cache.radius_and_derivative(params),
        cache.shape(params),
        cache.rho_z_grid(params, n_points),
        cache.rho_z_grid_unchecked(params, n_points),
        cache.neck(params),
        cache.geometry(params),
        cache.moller_coordinates(params),
        cache.scission_sigma1(params),
    ]


def test_invalid_coordinates_is_a_status() -> None:
    for res in _all_one_shot(BAD_K):
        assert res.status == tp.Status.invalid_coordinates


def test_wrong_count() -> None:
    with tp.Cache(tp.theta_grid(8)) as cache:
        assert cache.radius_grid(SYM_NECKED[:4]).status == tp.Status.wrong_param_count
        assert cache.radius_grid(SYM_NECKED + [0.0]).status == tp.Status.wrong_param_count
    assert tp.shape(SYM_NECKED + [0.0]).status == tp.Status.too_many_params
    assert tp.shape(SYM_NECKED[:4]).status == tp.Status.wrong_param_count
    assert tp.shape([]).status == tp.Status.wrong_param_count


def test_cache_rejections_raise() -> None:
    with pytest.raises(tp.TqsParamError) as err:
        tp.Cache(thetas=[])
    assert err.value.status == tp.Status.invalid_grid
    with pytest.raises(tp.TqsParamError) as err:
        tp.Cache(tp.theta_grid(8), max_params=4)
    assert err.value.status == tp.Status.invalid_init
    with pytest.raises(tp.TqsParamError) as err:
        tp.Cache(tp.theta_grid(8), max_params=6)
    assert err.value.status == tp.Status.too_many_params
    with pytest.raises(tp.TqsParamError) as err:
        tp.Cache([0.5, 1.5 * np.pi])
    assert err.value.status == tp.Status.invalid_grid


def test_closed_cache_raises() -> None:
    cache = tp.Cache(tp.theta_grid(8))
    cache.close()
    cache.close()
    with pytest.raises(tp.TqsParamError):
        cache.radius_grid(SYM_NECKED)
    with pytest.raises(tp.TqsParamError):
        cache.geometry(SYM_NECKED)


def test_non_1d_input_raises() -> None:
    with pytest.raises(tp.TqsParamError):
        tp.shape(np.zeros((5, 1)))


def test_inputs_are_converted_to_contiguous_float64() -> None:
    thetas = tp.theta_grid(16)
    base = [2.0, -1.0, 1.0, 0.0, 0.0]     # exact in float32 and as ints
    ref = tp.radius_grid(base, thetas)
    assert ref.ok
    variants = [
        tp.radius_grid([2, -1, 1, 0, 0], thetas),
        tp.radius_grid(np.array(base, dtype=np.float32), thetas),
        tp.radius_grid(np.array([2.0, 9.0, -1.0, 9.0, 1.0, 9.0, 0.0, 9.0, 0.0])[::2], thetas),
        tp.radius_grid(tuple(base), list(thetas)),
    ]
    for res in variants:
        assert res.ok
        assert _same_bits(res.radii, ref.radii)


def _numbers(res) -> list:
    """Every numeric output of a result object, flattened."""
    out = []
    for name, value in vars(res).items():
        if name == "status":
            continue
        if isinstance(value, np.ndarray):
            out.extend(value.tolist())
        elif isinstance(value, tuple):
            out.extend(float(v) for v in value)
        else:
            out.append(float(value))
    return out


@pytest.mark.parametrize("params", [BAD_K, NO_GEOM])
def test_failures_zero_fill_every_output(params) -> None:
    for res in _all_one_shot(params):
        if res.ok:
            continue
        assert all(v == 0.0 for v in _numbers(res)), type(res).__name__


def test_profile_floor() -> None:
    assert tp.N_POINTS_FLOOR == 3
    prof = tp.rho_z_grid(SYM_NECKED, 3)
    assert prof.ok and prof.rho[1] > 0.0
    assert tp.rho_z_grid(SYM_NECKED, 2).status == tp.Status.invalid_grid
    with tp.Cache(tp.theta_grid(8)) as cache:
        assert cache.rho_z_grid(SYM_NECKED, 2).status == tp.Status.invalid_grid


def _same_result(a, b) -> bool:
    if type(a) is not type(b) or a.status != b.status:
        return False
    for name, value in vars(a).items():
        other = getattr(b, name)
        if isinstance(value, np.ndarray):
            if not _same_bits(value, other):
                return False
        elif isinstance(value, tuple):
            if not _same_bits(np.array(value), np.array(other)):
                return False
        elif isinstance(value, float):
            if not _same_bits(np.array([value]), np.array([other])):
                return False
        elif value != other:
            return False
    return True


@pytest.mark.parametrize("params", [SYM_NECKED, PEAR, SEPARATED, BAD_K])
def test_one_shot_equals_cached_bits(params) -> None:
    thetas = tp.theta_grid(8)
    with tp.Cache(thetas) as cache:
        cached = _all_cached(cache, params)
    for a, b in zip(_all_one_shot(params, thetas), cached):
        assert _same_result(a, b), type(a).__name__


def test_at_thetas_equals_grid_bits() -> None:
    thetas = tp.theta_grid(32)
    with tp.Cache(thetas) as cache:
        for params in (SYM_NECKED, PEAR):
            grid = cache.radius_and_derivative(params)
            at = cache.radius_and_derivative_at_thetas(params, thetas)
            assert grid.ok and at.ok
            assert _same_bits(grid.radii, at.radii) and _same_bits(grid.dr_dtheta, at.dr_dtheta)


def test_cache_copies_its_thetas() -> None:
    thetas = tp.theta_grid(16)
    with tp.Cache(thetas) as cache:
        before = cache.radius_grid(ASYM_NECKED)
        thetas[:] = 0.3
        after = cache.radius_grid(ASYM_NECKED)
    assert before.ok and after.ok and _same_bits(before.radii, after.radii)


def test_cache_is_shared_between_threads() -> None:
    """Contract family 3: many threads computing on ONE cache."""
    n_threads, n_steps = 8, 100

    def shape_for(t: int, s: int) -> list:
        return [1.2 + 0.012 * s + 0.01 * t, -0.6 - 0.002 * s, 1.0, 0.05 * t - 0.2, 0.1 - 0.004 * s]

    def walk(cache: "tp.Cache", tid: int) -> list:
        out = []
        for step in range(n_steps):
            p = shape_for(tid, step)
            out.append((cache.radius_and_derivative(p), cache.rho_z_grid(p, 65), cache.geometry(p)))
        return out

    with tp.Cache(tp.theta_grid(48)) as cache:
        serial = [walk(cache, tid) for tid in range(n_threads)]
        threaded: list = [None] * n_threads

        def worker(tid: int) -> None:
            threaded[tid] = walk(cache, tid)

        pool = [threading.Thread(target=worker, args=(tid,)) for tid in range(n_threads)]
        for th in pool:
            th.start()
        for th in pool:
            th.join()

    n_valid = 0
    for tid in range(n_threads):
        assert threaded[tid] is not None
        for a, b in zip(threaded[tid], serial[tid]):
            assert all(_same_result(x, y) for x, y in zip(a, b))
            n_valid += int(all(r.ok for r in b))
    # Zero-filled rejections compare equal too: the test must not pass vacuously.
    assert n_valid > 0


def test_status_message_and_enum_coverage() -> None:
    for status in tp.Status:
        msg = tp.status_message(status)
        assert isinstance(msg, str) and msg and msg != "unknown status code"
    assert tp.status_message(999) == "unknown status code"
    assert tp.N_PARAMS == 5
    assert [int(k) for k in tp.MidKind] == [0, 1, 2, 3, 4, 5]
