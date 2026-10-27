# tqs-parameterization

Three-quadratic-surface (3QS) nuclear shape parameterization: a Fortran 2018
library with a C API and Python bindings. It takes Nix's five coordinates
`[sigma1, sigma2, sigma3, alpha2, alpha3]` and returns the axial surface:
R(theta) and dR/dtheta at caller thetas, the rho(z) profile, the poles, the
neck, the resolved geometry, Möller's grid coordinates as diagnostics, and the
scission locus. Everything is closed form: no Newton iteration, no
quadrature, no root search.

The library follows the two-tier shape parameterization contract shared with
[beta-parameterization](https://github.com/AleksanderAugustyn/beta-parameterization)
and [fos-parameterization](https://github.com/AleksanderAugustyn/fos-parameterization):
one-shot functions, plus a read-only cache that is built once and shared across
threads. Unlike those two it is fixed-arity: every call takes exactly five
parameters.

## The shape

J. R. Nix, Nucl. Phys. A 130 (1969) 241, eq. 4: three quadratic surfaces of
revolution, joined so that rho^2 and its z-derivative are continuous.

    rho^2 = a_i^2 - k_i (z - l_i)^2,     k_i = a_i^2 / c_i^2

Segment 1 is the south end spheroid, segment 2 the north one, and segment 3
the middle body between the joins z1 <= z2. The ends have k1, k2 > 0; the
middle has k3 = sigma2 of either sign.

| Coordinate | Definition | Meaning | Sphere |
|---|---|---|---|
| `sigma1` | (l2 - l1)/u | separation of the end-spheroid centres | 0 |
| `sigma2` | a3^2/c3^2 | curvature of the middle body; < 0 is a neck | 1 |
| `sigma3` | (a1^2/c1^2 + a2^2/c2^2)/2 | mean squared axis ratio of the ends | 1 |
| `alpha2` | (a1^2 - a2^2)/u^2 | end transverse-radius (mass) asymmetry | 0 |
| `alpha3` | a1^2/c1^2 - a2^2/c2^2 | end axis-ratio asymmetry | 0 |

with u = [(a1^2 + a2^2)/2]^(1/2). The sphere is `[0, 1, 1, 0, 0]`, and
`[0, s2, 1, 0, 0]` for any s2 (the middle has zero length). The library places
the centre of mass at the origin and works in R0 = 1 units (volume 4 pi/3);
Nix's alpha1 = (l1 + l2)/(2u) is reported, not taken. It has no notion of A, Z
or N.

The middle body is classified by sigma2, by a3^2, and by whether its waist l3
lies strictly inside (z1, z2):

| `mid_kind` | Condition | Shape |
|---|---|---|
| spheroid | sigma2 > 0 | convex middle |
| cylinder | sigma2 = 0 | needs alpha2 = 0 |
| one_sheet | sigma2 < 0, a3^2 > 0 | necked if the waist is inside, else a monotone flare |
| cone | sigma2 < 0, a3^2 = 0 | pinched (scission) if the waist is inside, else a frustum |
| two_sheet | sigma2 < 0, a3^2 < 0 | separated if the waist is inside, else a connected pear |

Every row is a valid shape. Pinched and separated shapes are flagged in the
geometry and refused by the R(theta) outputs (codes 102 and 103); the profile,
neck, geometry and Möller outputs evaluate them. Rejected outright are invalid
coordinates (code 100: k1 or k2 <= 0, |alpha2| >= 2, sigma1 < 0) and
coordinates with no smooth three-surface geometry (code 101).

Two rules worth knowing:

- **Degenerate curvature.** An end whose curvature equals the middle's
  (|k_i - sigma2| <= 1e-4) has no join of its own; it is accepted only as the
  two surfaces coinciding (a3^2 within 1e-3 of a_i^2), with a zero-length
  middle. Consequently the sphere and spheroids are isolated points on the
  line sigma2 = sigma3, alpha3 = 0: any sigma1 > 0 there is code 101.
- **The root rule is a guard.** The joins reduce to a quadratic in l3. If
  both roots gave a valid geometry the library would take the one nearer the
  end-centre midpoint; over the 110,565-point test box this never happens.

## Two tiers

1. **One-shot.** One call computes a shape and keeps nothing.
2. **Read-only cache.** Build a cache once on a primary theta set, then
   request any number of shapes against it. Nothing derived from a shape
   parameter is stored, so no call influences a later one, and any number of
   threads may compute on one cache at once.

The two tiers return the same status and bitwise-identical outputs.

| Output | Cached | One-shot |
|---|---|---|
| R(theta) at the primary thetas | `cache_radius_grid_s` | `compute_radius_grid_standalone_s` |
| R(theta) and dR/dtheta | `cache_radius_and_derivative_s` | `compute_radius_and_derivative_standalone_s` |
| R(theta) and dR/dtheta at caller thetas | `cache_radius_and_derivative_at_thetas_s` | (the one-shot forms take caller thetas) |
| pole radii, with the R(theta) verdict | `cache_shape_s` | `compute_shape_standalone_s` |
| rho(z) profile | `cache_rho_z_grid_s` | `compute_rho_z_grid_standalone_s` |
| rho(z) profile, separated shapes included | `cache_rho_z_grid_unchecked_s` | `compute_rho_z_grid_unchecked_standalone_s` |
| neck position and radius | `cache_neck_s` | `compute_neck_standalone_s` |
| resolved geometry and flags | `cache_geometry_s` | `compute_geometry_standalone_s` |
| Möller's q2, d, alpha_g, eps_f1, eps_f2 | `cache_moller_coordinates_s` | `compute_moller_coordinates_standalone_s` |
| scission locus sigma1_t, sigma1_star | `cache_scission_sigma1_s` | `compute_scission_sigma1_standalone_s` |

The C functions are `tqs_param_cache_<output>` and `tqs_param_<output>`; the
Python names are the methods of `Cache` and the module-level functions.
`q2` is the dimensionless <r^2 P2> over the volume; Möller's charge quadrupole
moment (his eq. 14) is Q2 = 2 Z R0^2 q2.

### Fortran

```fortran
use tqs_parameterization_mod, only: cache_t, cache_init_s, cache_radius_grid_s, &
        compute_geometry_standalone_s, tqs_geometry_t

type(cache_t) :: cache
type(tqs_geometry_t) :: geometry
real(kind = rk) :: thetas(64), radii(64)
integer(kind = ik) :: status

! One-shot
call compute_geometry_standalone_s([1.6_rk, -0.6_rk, 1.0_rk, 0.0_rk, 0.0_rk], geometry, status)

! Cached
call cache_init_s(cache, 5_ik, thetas, status)
call cache_radius_grid_s(cache, [1.6_rk, -0.6_rk, 1.0_rk, 0.3_rk, 0.2_rk], radii, status)
```

### C

```c
#include "tqs_parameterization.h"

const double params[TQS_PARAM_N_PARAMS] = {1.6, -0.6, 1.0, 0.3, 0.2};
int32_t status;
tqs_param_cache_t *cache = tqs_param_cache_create(TQS_PARAM_N_PARAMS, thetas, n_thetas, &status);
status = tqs_param_cache_radius_grid(cache, params, TQS_PARAM_N_PARAMS, radii, n_thetas);
tqs_param_geometry_t geometry;
tqs_param_geometry(params, TQS_PARAM_N_PARAMS, &geometry, &status);   /* one-shot */
tqs_param_cache_destroy(cache);
```

### Python

```python
import tqs_parameterization as tp

params = [1.6, -0.6, 1.0, 0.3, 0.2]
print(tp.geometry(params))                     # one-shot

with tp.Cache(tp.theta_grid(64)) as cache:     # cached
    res = cache.radius_and_derivative(params)
    if res.ok:
        print(res.radii, res.dr_dtheta)
```

## Status codes

Every fallible routine reports an integer status; `0` is success. On any
nonzero status every output is zero-filled.

| Code | Name | Meaning |
|---|---|---|
| 1 | `too_many_params` | `max_params > 5`, or a one-shot call with more than 5 parameters |
| 2 | `cache_not_initialized` | compute on a cache that was not built |
| 3 | `invalid_grid` | no thetas, a theta outside `[0, pi]`, a profile under 3 nodes, or an output array of the wrong size |
| 4 | `wrong_param_count` | not exactly 5 parameters |
| 5 | `invalid_init` | `max_params < 5` |
| 100 | `invalid_coordinates` | k1 or k2 <= 0, \|alpha2\| >= 2, or sigma1 < 0 |
| 101 | `no_geometry` | no smooth three-surface geometry for these coordinates |
| 102 | `not_star_convex` | R(theta) outputs only: pinched, or star-convexity margin below 1e-2 |
| 103 | `separated` | R(theta) outputs and the checked profile: two fragments |

Codes 1-5 are shared with the other shape libraries. Inputs must be finite:
the release build cannot detect NaN.

## Build and install

**Python.** 0.1.0 is not on PyPI yet. Build the wheel from a checkout:

    python -m pip install build
    python -m build --wheel --outdir dist_raw
    python -m pip install dist_raw/tqs_parameterization-*.whl

**CMake.** Fetch the library into a Fortran or C/C++ project:

```cmake
include(FetchContent)
FetchContent_Declare(
        tqs-parameterization
        GIT_REPOSITORY https://github.com/AleksanderAugustyn/tqs-parameterization.git
        GIT_TAG 0.1.0
)
set(TQS_PARAM_BUILD_TESTS OFF)
FetchContent_MakeAvailable(tqs-parameterization)

target_link_libraries(my_target PRIVATE TqsParameterization::tqs_parameterization)
```

Targets: `TqsParameterization::tqs_parameterization` (static, Fortran
modules), `TqsParameterization::tqs_parameterization_shared` (shared library)
and `TqsParameterization::tqs_parameterization_cxx` (the C header plus the
shared library).

**From a checkout.**

    cmake -S . -B build/release -DCMAKE_BUILD_TYPE=Release
    cmake --build build/release -j
    ctest --test-dir build/release --output-on-failure
    TQS_PARAM_LIB=$PWD/build/release/libtqs_parameterization.so PYTHONPATH=$PWD/python \
        python -m pytest python/tests -q

## Version pins

| Dependency | Version |
|---|---|
| [gcc-compiler-options](https://github.com/AleksanderAugustyn/gcc-compiler-options) | 2.0.0 |
| [fortran-foundations](https://github.com/AleksanderAugustyn/fortran-foundations) | 3.0.0 |

Both are fetched by CMake. The library needs gfortran 10 or newer; it is
developed with gfortran 13. Python 3.9 or newer, NumPy 1.21 or newer.

## Changelog and license

See the
[changelog](https://github.com/AleksanderAugustyn/tqs-parameterization/blob/master/CHANGELOG.md).
MIT license.
