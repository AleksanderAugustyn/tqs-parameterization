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

    !---------------------------------------------------------------------------
    ! Geometry record
    !---------------------------------------------------------------------------
    !> A resolved shape. Segment i is the surface rho^2 = a_sq(i) - k(i) (z - l(i))^2:
    !! 1 the south end spheroid on [z_south, z_join(1)], 3 the middle body on
    !! [z_join(1), z_join(2)], 2 the north end spheroid on [z_join(2), z_north].
    !! c(i) = sqrt(a_sq(i)/k(i)) are the end semi-axes along z. Lengths are in
    !! units of u until the resolve rescales them; u is Nix's scale.
    type, public :: geometry_t
        real(kind = rk) :: u = 0.0_rk
        real(kind = rk) :: a_sq(3) = 0.0_rk
        real(kind = rk) :: k(3) = 0.0_rk
        real(kind = rk) :: c(2) = 0.0_rk
        real(kind = rk) :: l(3) = 0.0_rk
        real(kind = rk) :: z_join(2) = 0.0_rk
        real(kind = rk) :: z_south = 0.0_rk
        real(kind = rk) :: z_north = 0.0_rk
        real(kind = rk) :: w = 0.0_rk            !! Half-width of the void, 0 unless is_separated
        integer(kind = ik) :: mid_kind = TQS_MID_NONE
        logical :: waist_inside = .false.
        logical :: is_separated = .false.
        logical :: is_pinched = .false.
    end type geometry_t

    public :: segment_integrals_s
    public :: body_integrals_s
    public :: solve_quadratic_s

contains

    !===========================================================================
    ! INTEGRALS
    !===========================================================================

    !> Moments of one surface rho^2 = a_sq - k (z - l)^2 over [z_lo, z_hi].
    !!
    !! Every integral is written as (z_hi - z_lo) * P(t_lo, t_hi) with
    !! t = z - l, never as a difference of primitives F(t_hi) - F(t_lo): the
    !! difference cancels catastrophically when |t| is large (sigma1 = 1e8
    !! loses every digit), the factored form does not. The moments in z follow
    !! from the moments in t by the binomial shift, which only multiplies.
    !!
    !! @param[in]  a_sq, k, l  The surface
    !! @param[in]  z_lo, z_hi  The interval; z_hi < z_lo gives the negated integrals
    !! @param[out] i0          int rho^2 dz
    !! @param[out] i1          int z rho^2 dz
    !! @param[out] i2          int z^2 rho^2 dz
    !! @param[out] i4          int rho^4 dz
    pure subroutine segment_integrals_s(a_sq, k, l, z_lo, z_hi, i0, i1, i2, i4)

        real(kind = rk), intent(in) :: a_sq, k, l, z_lo, z_hi
        real(kind = rk), intent(out) :: i0, i1, i2, i4

        real(kind = rk) :: d, ta, tb, s2, s4, j0, j1, j2

        d = z_hi - z_lo
        ta = z_lo - l
        tb = z_hi - l
        s2 = tb * tb + tb * ta + ta * ta
        s4 = tb**4 + tb**3 * ta + tb * tb * ta * ta + tb * ta**3 + ta**4

        j0 = d * (a_sq - k * s2 / 3.0_rk)
        j1 = d * (0.5_rk * a_sq * (tb + ta) - 0.25_rk * k * (tb + ta) * (tb * tb + ta * ta))
        j2 = d * (a_sq * s2 / 3.0_rk - k * s4 / 5.0_rk)

        i0 = j0
        i1 = j1 + l * j0
        i2 = j2 + 2.0_rk * l * j1 + l * l * j0
        i4 = d * (a_sq * a_sq - 2.0_rk * a_sq * k * s2 / 3.0_rk + k * k * s4 / 5.0_rk)

    end subroutine segment_integrals_s

    !> Moments of the whole body: end 1, the middle (minus the void when
    !! separated), end 2, summed in that order.
    !!
    !! @param[in]  geom            Resolved geometry, any units
    !! @param[out] i0, i1, i2, i4  As segment_integrals_s, over the body
    pure subroutine body_integrals_s(geom, i0, i1, i2, i4)

        type(geometry_t), intent(in) :: geom
        real(kind = rk), intent(out) :: i0, i1, i2, i4

        real(kind = rk) :: p0, p1, p2, p4

        call segment_integrals_s(geom%a_sq(1), geom%k(1), geom%l(1), geom%z_south, &
                geom%z_join(1), i0, i1, i2, i4)

        if (geom%is_separated) then
            call segment_integrals_s(geom%a_sq(3), geom%k(3), geom%l(3), geom%z_join(1), &
                    geom%l(3) - geom%w, p0, p1, p2, p4)
            call accumulate_s(p0, p1, p2, p4, i0, i1, i2, i4)
            call segment_integrals_s(geom%a_sq(3), geom%k(3), geom%l(3), &
                    geom%l(3) + geom%w, geom%z_join(2), p0, p1, p2, p4)
        else
            call segment_integrals_s(geom%a_sq(3), geom%k(3), geom%l(3), geom%z_join(1), &
                    geom%z_join(2), p0, p1, p2, p4)
        end if
        call accumulate_s(p0, p1, p2, p4, i0, i1, i2, i4)

        call segment_integrals_s(geom%a_sq(2), geom%k(2), geom%l(2), geom%z_join(2), &
                geom%z_north, p0, p1, p2, p4)
        call accumulate_s(p0, p1, p2, p4, i0, i1, i2, i4)

    end subroutine body_integrals_s

    pure subroutine accumulate_s(p0, p1, p2, p4, i0, i1, i2, i4)
        real(kind = rk), intent(in) :: p0, p1, p2, p4
        real(kind = rk), intent(inout) :: i0, i1, i2, i4
        i0 = i0 + p0
        i1 = i1 + p1
        i2 = i2 + p2
        i4 = i4 + p4
    end subroutine accumulate_s

    !===========================================================================
    ! QUADRATIC
    !===========================================================================

    !> Real roots of a x^2 + b x + c = 0 in the stable form.
    !!
    !! q = -(b + sign(b) sqrt(disc))/2, roots q/a and c/q, so neither root is
    !! a difference of nearly equal numbers. A discriminant negative by less
    !! than QUADRATIC_DISC_TOL * max(b^2, |4 a c|) is clamped to zero: on the
    !! coincident family the exact discriminant is zero for every parameter
    !! value and roundoff gives either sign.
    !!
    !! @param[in]  a, b, c  Coefficients
    !! @param[out] n_roots  2; 1 (double root, or linear when a = 0); 0 (no real
    !!                      root, or a = b = 0 with c /= 0); -1 (a = b = c = 0,
    !!                      identically satisfied)
    !! @param[out] r1, r2   The roots; r2 = r1 when n_roots = 1, both 0 when n_roots <= 0
    pure subroutine solve_quadratic_s(a, b, c, n_roots, r1, r2)

        real(kind = rk), intent(in) :: a, b, c
        integer(kind = ik), intent(out) :: n_roots
        real(kind = rk), intent(out) :: r1, r2

        real(kind = rk) :: disc, scale, q

        r1 = 0.0_rk
        r2 = 0.0_rk

        if (.not. (abs(a) > 0.0_rk)) then
            if (abs(b) > 0.0_rk) then
                n_roots = 1_ik
                r1 = -c / b
                r2 = r1
            else if (abs(c) > 0.0_rk) then
                n_roots = 0_ik
            else
                n_roots = -1_ik
            end if
            return
        end if

        disc = b * b - 4.0_rk * a * c
        if (disc > 0.0_rk) then
            q = -0.5_rk * (b + sign(sqrt(disc), b))
            n_roots = 2_ik
            r1 = q / a
            r2 = c / q
            return
        end if

        scale = max(b * b, abs(4.0_rk * a * c))
        if (-disc <= QUADRATIC_DISC_TOL * scale) then
            n_roots = 1_ik
            r1 = -0.5_rk * b / a
            r2 = r1
        else
            n_roots = 0_ik
        end if

    end subroutine solve_quadratic_s

end module tqs_workers_mod
