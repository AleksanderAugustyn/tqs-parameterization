# tqs-parameterization

A Fortran library for the three-quadratic-surface (3QS) family of axially
symmetric nuclear shapes of Nix (Nucl. Phys. A 130 (1969) 241), in Nix's own
coordinates `[sigma1, sigma2, sigma3, alpha2, alpha3]`. The inversion to
geometry is closed form; the library emits R(theta), dR/dtheta, the rho(z)
profile, poles, neck, the resolved geometry, Moller's grid coordinates as
diagnostics and the scission locus, through Fortran, C and Python
interfaces on the two-tier shape-parameterization contract.
