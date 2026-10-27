!> Named fixtures of the tqs test suites and the independent quadrature
!! oracle. Values are [sigma1, sigma2, sigma3, alpha2, alpha3]; the
!! classification and cross-check numbers live in the plan's fixture table.
module tqs_test_fixtures_mod

    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use mathematical_utilities_mod, only: compute_gauss_legendre_quadrature_s
    use tqs_workers_mod, only: geometry_t, segment_rho_s
    use tqs_resolve_mod, only: resolve_s
    use test_utils_mod, only: assert_int_eq

    implicit none

    private

    public :: SPHERE, SPHERE_B, SPHEROID_K2, SPHEROID_S2, SYM_NECKED, INITIAL, &
            ASYM_NECKED, SPHEROID_MID, CYLINDER, PEAR, PINCHED, SEPARATED, &
            SEPARATED_ASYM, STAR_FAIL, DEGENERATE_ONE, BAD_K, BAD_A2, BAD_S1, &
            NO_GEOM, NO_CYL
    public :: pv, resolve_ok, gl_body_moments

    real(kind = rk), parameter :: SPHERE(5) = [0.0_rk, 1.0_rk, 1.0_rk, 0.0_rk, 0.0_rk]
    real(kind = rk), parameter :: SPHERE_B(5) = [0.0_rk, -1.0_rk, 1.0_rk, 0.0_rk, 0.0_rk]
    real(kind = rk), parameter :: SPHEROID_K2(5) = [0.0_rk, 2.0_rk, 2.0_rk, 0.0_rk, 0.0_rk]
    real(kind = rk), parameter :: SPHEROID_S2(5) = [0.0_rk, -0.5_rk, 2.0_rk, 0.0_rk, 0.0_rk]
    real(kind = rk), parameter :: SYM_NECKED(5) = [1.6_rk, -0.6_rk, 1.0_rk, 0.0_rk, 0.0_rk]
    real(kind = rk), parameter :: INITIAL(5) = [2.4_rk, -1.0_rk, 1.0_rk, 0.0_rk, 0.0_rk]
    real(kind = rk), parameter :: ASYM_NECKED(5) = [1.6_rk, -0.6_rk, 1.0_rk, 0.3_rk, 0.2_rk]
    real(kind = rk), parameter :: SPHEROID_MID(5) = [0.8_rk, 0.5_rk, 1.0_rk, 0.2_rk, 0.1_rk]
    real(kind = rk), parameter :: CYLINDER(5) = [0.6_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.3_rk]
    real(kind = rk), parameter :: PEAR(5) = [1.0_rk, -1.0_rk, 1.0_rk, -1.4_rk, 0.0_rk]
    real(kind = rk), parameter :: PINCHED(5) = [2.0_rk * sqrt(2.0_rk), -1.0_rk, 1.0_rk, 0.0_rk, 0.0_rk]
    real(kind = rk), parameter :: SEPARATED(5) = [3.0_rk, -1.0_rk, 1.0_rk, 0.0_rk, 0.0_rk]
    real(kind = rk), parameter :: SEPARATED_ASYM(5) = [3.4_rk, -0.6_rk, 1.0_rk, 0.3_rk, 0.2_rk]
    real(kind = rk), parameter :: STAR_FAIL(5) = [3.0_rk, -0.6_rk, 1.0_rk, 0.9_rk, 0.0_rk]
    real(kind = rk), parameter :: DEGENERATE_ONE(5) = [0.5_rk, 2.0_rk, 1.5_rk, -0.5_rk, 1.0_rk]
    real(kind = rk), parameter :: BAD_K(5) = [0.0_rk, 1.0_rk, 0.5_rk, 0.0_rk, 1.0_rk]
    real(kind = rk), parameter :: BAD_A2(5) = [0.5_rk, -0.6_rk, 1.0_rk, 2.0_rk, 0.0_rk]
    real(kind = rk), parameter :: BAD_S1(5) = [-0.1_rk, -0.6_rk, 1.0_rk, 0.0_rk, 0.0_rk]
    real(kind = rk), parameter :: NO_GEOM(5) = [0.3_rk, 1.5_rk, 1.5_rk, 0.0_rk, 0.0_rk]
    real(kind = rk), parameter :: NO_CYL(5) = [0.3_rk, 0.0_rk, 1.0_rk, 0.2_rk, 0.0_rk]

    integer(kind = ik), parameter :: N_GL = 64_ik

contains

    !> A parameter vector built element by element (a constructor with
    !! variable elements is an array temporary, an error in Debug).
    pure function pv(s1, s2, s3, a2, a3) result(p)
        real(kind = rk), intent(in) :: s1, s2, s3, a2, a3
        real(kind = rk) :: p(5)
        p(1) = s1
        p(2) = s2
        p(3) = s3
        p(4) = a2
        p(5) = a3
    end function pv

    !> resolve_s with the status asserted to be 0.
    subroutine resolve_ok(params, geom, label)
        real(kind = rk), intent(in) :: params(:)
        type(geometry_t), intent(out) :: geom
        character(len = *), intent(in) :: label
        integer(kind = ik) :: status
        call resolve_s(params, geom, status)
        call assert_int_eq(status, 0_ik, label // ': resolve status 0')
    end subroutine resolve_ok

    !> Independent oracle: 64-node Gauss-Legendre per segment of the profile
    !! from segment_rho_s, the void excluded when separated. Returns
    !! vol = pi int rho^2, z_moment = int z rho^2, q2 = 3/4 int (z^2 rho^2 - rho^4/4).
    subroutine gl_body_moments(geom, vol, z_moment, q2)
        type(geometry_t), intent(in) :: geom
        real(kind = rk), intent(out) :: vol, z_moment, q2
        real(kind = rk) :: x(N_GL), wt(N_GL), m0, m1, mq
        call compute_gauss_legendre_quadrature_s(N_GL, x, wt)
        m0 = 0.0_rk
        m1 = 0.0_rk
        mq = 0.0_rk
        call add_segment(1_ik, geom%z_south, geom%z_join(1))
        if (geom%is_separated) then
            call add_segment(3_ik, geom%z_join(1), geom%l(3) - geom%w)
            call add_segment(3_ik, geom%l(3) + geom%w, geom%z_join(2))
        else
            call add_segment(3_ik, geom%z_join(1), geom%z_join(2))
        end if
        call add_segment(2_ik, geom%z_join(2), geom%z_north)
        vol = PI_C * m0
        z_moment = m1
        q2 = 0.75_rk * mq
    contains
        subroutine add_segment(seg, z_lo, z_hi)
            integer(kind = ik), intent(in) :: seg
            real(kind = rk), intent(in) :: z_lo, z_hi
            real(kind = rk) :: half, mid, z, rho, drho, rho_sq, wz
            integer(kind = ik) :: j
            half = 0.5_rk * (z_hi - z_lo)
            mid = 0.5_rk * (z_hi + z_lo)
            do j = 1_ik, N_GL
                z = mid + half * x(j)
                wz = half * wt(j)
                call segment_rho_s(geom%a_sq(seg), geom%k(seg), geom%l(seg), z, rho, drho)
                rho_sq = rho * rho
                m0 = m0 + wz * rho_sq
                m1 = m1 + wz * z * rho_sq
                mq = mq + wz * (z * z * rho_sq - 0.25_rk * rho_sq * rho_sq)
            end do
        end subroutine add_segment
    end subroutine gl_body_moments

end module tqs_test_fixtures_mod
