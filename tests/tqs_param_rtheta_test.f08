!> R(theta) by closed-form ray intersection, dR/dtheta by implicit
!! differentiation, and the closed-form star-convexity margin against a
!! dense sampled minimum.
program tqs_param_rtheta_test
    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use tqs_workers_mod, only: geometry_t, rho_at_z_s, ray_radius_s, star_margin_f
    use tqs_test_fixtures_mod, only: SPHERE, SPHERE_B, SPHEROID_K2, SPHEROID_S2, &
            SYM_NECKED, INITIAL, ASYM_NECKED, SPHEROID_MID, CYLINDER, PEAR, PINCHED, &
            STAR_FAIL, DEGENERATE_ONE, pv, resolve_ok
    use test_utils_mod, only: assert_true, assert_abs_close, test_summary
    implicit none

    integer(kind = ik), parameter :: N_STAR = 8_ik
    real(kind = rk), parameter :: STAR_SET(5, N_STAR) = reshape([SYM_NECKED, INITIAL, &
            ASYM_NECKED, SPHEROID_MID, CYLINDER, PEAR, SPHERE, SPHEROID_S2], [5, N_STAR])
    character(len = 12), parameter :: STAR_NAMES(N_STAR) = [character(len = 12) :: &
            'SYM_NECKED', 'INITIAL', 'ASYM_NECKED', 'SPHEROID_MID', 'CYLINDER', 'PEAR', &
            'SPHERE', 'SPHEROID_S2']

    call test_t1_profile_points()
    call test_t2_derivative()
    call test_t3_asymptotic_cone()
    call test_t4_poles()
    call test_t5_star_gate()
    call test_t6_margin_cross_checks()

    call test_summary()

contains

    !> R(theta) at an angle, through ray_radius_s.
    subroutine radius_at(g, theta, r, dr, found)
        type(geometry_t), intent(in) :: g
        real(kind = rk), intent(in) :: theta
        real(kind = rk), intent(out) :: r, dr
        logical, intent(out) :: found
        call ray_radius_s(g, cos(theta), sin(theta), r, dr, found)
    end subroutine radius_at

    subroutine test_t1_profile_points()
        integer(kind = ik), parameter :: N = 401_ik
        type(geometry_t) :: g
        real(kind = rk) :: z, rho, drho, r, dr, worst
        logical :: found, all_found
        integer(kind = ik) :: f, j
        do f = 1_ik, N_STAR
            call resolve_ok(STAR_SET(:, f), g, 'T1 ' // trim(STAR_NAMES(f)))
            worst = 0.0_rk
            all_found = .true.
            do j = 1_ik, N
                z = g%z_south + (g%z_north - g%z_south) * (real(j, rk) / real(N + 1_ik, rk))
                call rho_at_z_s(g, z, rho, drho)
                if (.not. (rho > 0.0_rk)) cycle
                call radius_at(g, atan2(rho, z), r, dr, found)
                all_found = all_found .and. found
                worst = max(worst, abs(r - hypot(z, rho)))
            end do
            call assert_true(all_found, 'T1 ' // trim(STAR_NAMES(f)) // ': every ray hits')
            call assert_true(worst <= 1.0e-13_rk, 'T1 ' // trim(STAR_NAMES(f)) // ': R = hypot(z, rho)')
        end do
    end subroutine test_t1_profile_points

    !> dR/dtheta against central differences at one angle; returns .true. if it passes.
    logical function derivative_ok(g, theta) result(ok)
        type(geometry_t), intent(in) :: g
        real(kind = rk), intent(in) :: theta
        real(kind = rk), parameter :: H = 1.0e-6_rk
        real(kind = rk) :: r, dr, rp, rm, dummy, fd
        logical :: f0, fp, fm
        call radius_at(g, theta, r, dr, f0)
        call radius_at(g, theta + H, rp, dummy, fp)
        call radius_at(g, theta - H, rm, dummy, fm)
        fd = (rp - rm) / (2.0_rk * H)
        ok = f0 .and. fp .and. fm .and. abs(dr - fd) <= 1.0e-6_rk * (1.0_rk + abs(dr))
    end function derivative_ok

    subroutine derivative_sweep(params, label)
        real(kind = rk), intent(in) :: params(5)
        character(len = *), intent(in) :: label
        type(geometry_t) :: g
        real(kind = rk) :: rho, drho, th_join
        integer(kind = ik) :: i, n_bad, side
        call resolve_ok(params, g, label)
        n_bad = 0_ik
        do i = 1_ik, 181_ik
            if (.not. derivative_ok(g, real(i, rk) * PI_C / 182.0_rk)) n_bad = n_bad + 1_ik
        end do
        do side = 1_ik, 2_ik
            call rho_at_z_s(g, g%z_join(side), rho, drho)
            th_join = atan2(rho, g%z_join(side))
            if (.not. derivative_ok(g, th_join - 1.0e-4_rk)) n_bad = n_bad + 1_ik
            if (.not. derivative_ok(g, th_join + 1.0e-4_rk)) n_bad = n_bad + 1_ik
        end do
        call assert_true(n_bad == 0_ik, label // ': dR/dtheta matches central differences')
    end subroutine derivative_sweep

    subroutine test_t2_derivative()
        call derivative_sweep(ASYM_NECKED, 'T2 ASYM_NECKED')
        call derivative_sweep(PEAR, 'T2 PEAR')
    end subroutine test_t2_derivative

    subroutine test_t3_asymptotic_cone()
        type(geometry_t) :: g
        real(kind = rk) :: th0, th, r, dr, rho, drho
        logical :: found
        integer(kind = ik) :: j
        character(len = 24) :: tag
        call resolve_ok(SYM_NECKED, g, 'T3 SYM_NECKED')
        th0 = atan(sqrt(0.6_rk))
        do j = -1_ik, 1_ik
            write(tag, '(A,I0)') 'T3 offset ', j
            th = th0 + real(j, rk) * 1.0e-3_rk
            call radius_at(g, th, r, dr, found)
            call assert_true(found, trim(tag) // ': found')
            call rho_at_z_s(g, r * cos(th), rho, drho)
            call assert_abs_close(rho, r * sin(th), 1.0e-12_rk, trim(tag) // ': the hit lies on the profile')
            call assert_true(derivative_ok(g, th), trim(tag) // ': dR/dtheta matches central differences')
        end do
    end subroutine test_t3_asymptotic_cone

    subroutine pole_check(params, label)
        real(kind = rk), intent(in) :: params(5)
        character(len = *), intent(in) :: label
        type(geometry_t) :: g
        real(kind = rk) :: r, dr
        logical :: found
        call resolve_ok(params, g, label)
        call ray_radius_s(g, 1.0_rk, 0.0_rk, r, dr, found)
        call assert_true(found, label // ': north pole found')
        call assert_abs_close(r, g%z_north, 1.0e-14_rk, label // ': R(0) = r_north')
        call assert_abs_close(dr, 0.0_rk, 1.0e-15_rk, label // ': dR/dtheta(0) = 0')
        call ray_radius_s(g, -1.0_rk, 0.0_rk, r, dr, found)
        call assert_true(found, label // ': south pole found')
        call assert_abs_close(r, -g%z_south, 1.0e-14_rk, label // ': R(pi) = r_south')
        call assert_abs_close(dr, 0.0_rk, 1.0e-15_rk, label // ': dR/dtheta(pi) = 0')
    end subroutine pole_check

    subroutine test_t4_poles()
        call pole_check(SPHERE, 'T4 SPHERE')
        call pole_check(SPHERE_B, 'T4 SPHERE_B')
        call pole_check(SPHEROID_K2, 'T4 SPHEROID_K2')
        call pole_check(SPHEROID_S2, 'T4 SPHEROID_S2')
        call pole_check(SYM_NECKED, 'T4 SYM_NECKED')
        call pole_check(INITIAL, 'T4 INITIAL')
        call pole_check(ASYM_NECKED, 'T4 ASYM_NECKED')
        call pole_check(SPHEROID_MID, 'T4 SPHEROID_MID')
        call pole_check(CYLINDER, 'T4 CYLINDER')
        call pole_check(PEAR, 'T4 PEAR')
        call pole_check(PINCHED, 'T4 PINCHED')
        call pole_check(STAR_FAIL, 'T4 STAR_FAIL')
        call pole_check(DEGENERATE_ONE, 'T4 DEGENERATE_ONE')
    end subroutine test_t4_poles

    !> The densely sampled minimum of rho - z drho/dz over interior nodes with rho > 0.
    real(kind = rk) function sampled_margin(g) result(m)
        type(geometry_t), intent(in) :: g
        integer(kind = ik), parameter :: N = 200001_ik
        real(kind = rk) :: z, rho, drho
        integer(kind = ik) :: j
        m = huge(1.0_rk)
        do j = 1_ik, N
            z = g%z_south + (g%z_north - g%z_south) * (real(j, rk) / real(N + 1_ik, rk))
            call rho_at_z_s(g, z, rho, drho)
            if (rho > 0.0_rk) m = min(m, rho - z * drho)
        end do
    end function sampled_margin

    subroutine margin_check(params, label, closed_out)
        real(kind = rk), intent(in) :: params(5)
        character(len = *), intent(in) :: label
        real(kind = rk), intent(out), optional :: closed_out
        type(geometry_t) :: g
        real(kind = rk) :: closed, sampled
        call resolve_ok(params, g, label)
        closed = star_margin_f(g)
        sampled = sampled_margin(g)
        call assert_true(closed <= sampled + 1.0e-12_rk, label // ': closed margin <= sampled')
        call assert_true(sampled - closed <= 1.0e-4_rk, label // ': sampled - closed <= 1e-4')
        if (present(closed_out)) closed_out = closed
    end subroutine margin_check

    subroutine test_t5_star_gate()
        real(kind = rk) :: m, p(5)
        integer(kind = ik) :: f, j
        character(len = 32) :: tag
        do f = 1_ik, N_STAR
            call margin_check(STAR_SET(:, f), 'T5 ' // trim(STAR_NAMES(f)))
        end do
        call margin_check(STAR_FAIL, 'T5 STAR_FAIL', m)
        call assert_abs_close(m, -0.5123_rk, 1.0e-4_rk, 'T5 STAR_FAIL margin')
        ! Pear family: waist-outside one_sheet and two_sheet middles.
        do j = 0_ik, 11_ik
            p = pv(1.0_rk, -1.0_rk, 1.0_rk, -1.45_rk + 0.05_rk * real(j, rk), 0.0_rk)
            write(tag, '(A,F6.2)') 'T5 pear alpha2 = ', p(4)
            call margin_check(p, trim(tag))
        end do
    end subroutine test_t5_star_gate

    subroutine test_t6_margin_cross_checks()
        type(geometry_t) :: g
        call resolve_ok(SYM_NECKED, g, 'T6 SYM_NECKED')
        call assert_abs_close(star_margin_f(g), 0.6278_rk, 1.0e-4_rk, 'T6 SYM_NECKED margin')
        call resolve_ok(ASYM_NECKED, g, 'T6 ASYM_NECKED')
        call assert_abs_close(star_margin_f(g), 0.5791_rk, 1.0e-4_rk, 'T6 ASYM_NECKED margin')
        call resolve_ok(PEAR, g, 'T6 PEAR')
        call assert_abs_close(star_margin_f(g), 0.9932_rk, 1.0e-4_rk, 'T6 PEAR margin')
        call resolve_ok(SPHERE, g, 'T6 SPHERE')
        call assert_abs_close(star_margin_f(g), 1.0_rk, 1.0e-4_rk, 'T6 SPHERE margin')
    end subroutine test_t6_margin_cross_checks

end program tqs_param_rtheta_test
