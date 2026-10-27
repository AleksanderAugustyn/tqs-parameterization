# Changelog

All notable changes to this project are documented here. Versions follow
the git tags.

## 0.1.0 — unreleased

First release: Nix's three-quadratic-surface shapes in Nix's coordinates
`[sigma1, sigma2, sigma3, alpha2, alpha3]`, on the two-tier
shape-parameterization contract as its first fixed-arity library.

- Closed-form resolve: end spheroids from the coordinates, the two smooth
  joins as a quadratic in the waist position l3 (stable form), volume scale
  and centre of mass from factored segment integrals. No iteration, no
  quadrature, no root search.
- Ten outputs in both tiers, bitwise identical between them: R(theta),
  R(theta) with dR/dtheta (at the cache's thetas and at caller thetas),
  pole radii, the checked and unchecked rho(z) profile, the neck, the
  resolved geometry with its middle-body kind and flags, Möller's q2, d,
  alpha_g, eps_f1, eps_f2, and the scission locus (Möller's sigma1_t and the
  a3^2 = 0 separation sigma1_star).
- Status codes 100 (invalid coordinates), 101 (no geometry), 102 (not
  star-convex or pinched, R(theta) outputs only), 103 (separated).
- C API (`tqs_param_`, `tqs_param_geometry_t`) and a ctypes Python package
  (`tqs_parameterization`).

Calls made in this release (owner-visible):

- An output array of the wrong extent is code 3, not a library code; a
  negative output size from C is the same.
- The rho(z) profile floor is 3 nodes (the poles and one interior node).
- `scission_sigma1` validates sigma1 like every other output, though its value
  does not depend on it.
- `geometry` is one record on every surface (Fortran `tqs_geometry_t`, C
  `tqs_param_geometry_t`, Python `Geometry`).
- A cylinder middle (sigma2 = 0) has l3 = (l1 + l2)/2.
- Degenerate curvature: |k_i - sigma2| <= 1e-4 counts as equal and is
  accepted only as surface coincidence, |a3^2 - a_i^2| <= 1e-3 max(1, a_i^2);
  both ends degenerate only on the exact spheroid line sigma1 = alpha2 = 0. A
  discriminant negative by under 1e-12 of its scale is a double root.
- Where the end spheroids are tangent the joins coincide; z1 above z2 by at
  most 1e-10 (roundoff) is one join, not a rejection.
- A cone middle (|a3^2| <= 1e-12) stores a3^2 = 0, so a pinch point is rho = 0
  exactly.
- Goldens live in the Python suite, captured through the C API.
- There is no standalone at-thetas form: the one-shot radius forms take caller
  thetas.
