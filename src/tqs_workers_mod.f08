!> Arithmetic on a resolved three-quadratic-surface (3QS) geometry.
!!
!! Owns the constants and enums of the library, the geometry record, and
!! every kernel that works on a resolved geometry: segment and body
!! integrals, the stable quadratic, the profile, the neck, ray intersection,
!! the star-convexity margin and the Moller coordinates. Everything is closed
!! form: no Newton iteration, no quadrature, no root search.
module tqs_workers_mod

    use precision_utilities_mod, only: ik, rk

    implicit none

    private

    !---------------------------------------------------------------------------
    ! Fixed arity and resolution floor
    !---------------------------------------------------------------------------
    !> params = [sigma1, sigma2, sigma3, alpha2, alpha3]; positions are permanent.
    integer(kind = ik), parameter, public :: TQS_N_PARAMS = 5_ik
    !> The two poles and one interior node: a sanity bound, the profile is closed form.
    integer(kind = ik), parameter, public :: TQS_N_POINTS_FLOOR = 3_ik

    !---------------------------------------------------------------------------
    ! Tolerances
    !---------------------------------------------------------------------------
    !> Minimum of rho - z drho/dz for the R(theta) outputs (representability, as fos-param).
    real(kind = rk), parameter, public :: STAR_CONVEXITY_MARGIN = 1.0e-2_rk
    !> |a3^2| in u^2 units at or below this counts as zero: roundoff at the
    !! scission locus must not flip the cone / two_sheet split.
    real(kind = rk), parameter, public :: TQS_PINCH_TOL = 1.0e-12_rk
    !> |k_i - sigma2| at or below this counts as equal (degenerate curvature).
    real(kind = rk), parameter, public :: TQS_DEGENERATE_TOL = 1.0e-4_rk
    !> |a3^2 - a_i^2| <= tol * max(1, a_i^2) is accepted as surface coincidence.
    real(kind = rk), parameter, public :: TQS_COINCIDENCE_TOL = 1.0e-3_rk
    !> A discriminant negative by less than this fraction of its scale is a
    !! double root: on the coincident family the exact discriminant is zero.
    real(kind = rk), parameter :: QUADRATIC_DISC_TOL = 1.0e-12_rk

    !---------------------------------------------------------------------------
    ! Middle-body kinds
    !---------------------------------------------------------------------------
    integer(kind = ik), parameter, public :: TQS_MID_NONE = 0_ik
    integer(kind = ik), parameter, public :: TQS_MID_SPHEROID = 1_ik
    integer(kind = ik), parameter, public :: TQS_MID_CYLINDER = 2_ik
    integer(kind = ik), parameter, public :: TQS_MID_ONE_SHEET = 3_ik
    integer(kind = ik), parameter, public :: TQS_MID_CONE = 4_ik
    integer(kind = ik), parameter, public :: TQS_MID_TWO_SHEET = 5_ik

    !---------------------------------------------------------------------------
    ! Library status codes (100+); the shared 0-5 live in shape_core_mod
    !---------------------------------------------------------------------------
    !! Owned here because the resolve module raises them and cannot depend on
    !! the front module; the front module re-exports them.
    integer(kind = ik), parameter, public :: TQS_ERROR_INVALID_COORDINATES = 100_ik
    integer(kind = ik), parameter, public :: TQS_ERROR_NO_GEOMETRY = 101_ik
    integer(kind = ik), parameter, public :: TQS_ERROR_NOT_STAR_CONVEX = 102_ik
    integer(kind = ik), parameter, public :: TQS_ERROR_SEPARATED = 103_ik

end module tqs_workers_mod
