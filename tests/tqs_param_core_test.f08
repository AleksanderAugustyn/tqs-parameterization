!> Core physics of the resolved shape: the sphere and spheroids by both
!! routes, C1 joins, volume and centre of mass against quadrature, the
!! analytic profile derivative, and the profile conventions at voids and
!! pinch points.
program tqs_param_core_test
    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use tqs_workers_mod, only: geometry_t, segment_rho_s, rho_at_z_s, profile_s, &
            poles_s, moller_coordinates_s
    use tqs_test_fixtures_mod, only: SPHERE, SPHERE_B, SPHEROID_K2, SPHEROID_S2, &
            SYM_NECKED, INITIAL, ASYM_NECKED, SPHEROID_MID, CYLINDER, PEAR, PINCHED, &
            SEPARATED, SEPARATED_ASYM, STAR_FAIL, DEGENERATE_ONE, resolve_ok, &
            gl_body_moments
    use test_utils_mod, only: assert_true, assert_abs_close, is_zero_f, test_summary
    implicit none

    call test_c1_sphere()
    call test_c2_spheroids()
    call test_c3_join_continuity()
    call test_c4_volume_and_com()
    call test_c5_profile_derivative()
    call test_c6_conventions()

    call test_summary()

contains

    subroutine sphere_checks(params, label, z, rho, drho)
        real(kind = rk), intent(in) :: params(5)
        character(len = *), intent(in) :: label
        real(kind = rk), intent(out) :: z(101), rho(101), drho(101)
        type(geometry_t) :: geom
        real(kind = rk) :: q2, d, ag, e1, e2, worst
        integer(kind = ik) :: j
        call resolve_ok(params, geom, label)
        call assert_abs_close(geom%u, 1.0_rk, 1.0e-14_rk, label // ': u = 1')
        call profile_s(geom, z, rho, drho)
        worst = 0.0_rk
        do j = 2_ik, 100_ik
            worst = max(worst, abs(rho(j)**2 - (1.0_rk - z(j)**2)))
        end do
        call assert_true(worst <= 1.0e-14_rk, label // ': rho^2 = 1 - z^2 at interior nodes')
        call assert_true(is_zero_f(rho(1)) .and. is_zero_f(drho(1)) .and. is_zero_f(rho(101)) &
                .and. is_zero_f(drho(101)), label // ': (0, 0) at both poles')
        call moller_coordinates_s(geom, q2, d, ag, e1, e2)
        call assert_abs_close(q2, 0.0_rk, 1.0e-14_rk, label // ': q2 = 0')
    end subroutine sphere_checks

    subroutine test_c1_sphere()
        real(kind = rk) :: za(101), ra(101), da(101), zb(101), rb(101), db(101)
        call sphere_checks(SPHERE, 'C1 SPHERE', za, ra, da)
        call sphere_checks(SPHERE_B, 'C1 SPHERE_B', zb, rb, db)
        call assert_true(maxval(abs(za - zb)) <= 1.0e-15_rk .and. maxval(abs(ra - rb)) <= 1.0e-15_rk &
                .and. maxval(abs(da - db)) <= 1.0e-15_rk, 'C1 the two sphere routes give one profile')
    end subroutine test_c1_sphere

    subroutine test_c2_spheroids()
        type(geometry_t) :: ga, gb
        real(kind = rk) :: za(101), ra(101), da(101), zb(101), rb(101), db(101), rn, rs
        call resolve_ok(SPHEROID_K2, ga, 'C2 SPHEROID_K2')
        call resolve_ok(SPHEROID_S2, gb, 'C2 SPHEROID_S2')
        call assert_abs_close(ga%a_sq(1) * ga%c(1), 1.0_rk, 1.0e-14_rk, 'C2 K2 a^2 c = 1')
        call assert_abs_close(gb%a_sq(1) * gb%c(1), 1.0_rk, 1.0e-14_rk, 'C2 S2 a^2 c = 1')
        call assert_abs_close(ga%k(1), 2.0_rk, 0.0_rk, 'C2 K2 k = 2')
        call assert_abs_close(gb%k(1), 2.0_rk, 0.0_rk, 'C2 S2 k = 2')
        call poles_s(ga, rn, rs)
        call assert_abs_close(rn, 0.793701_rk, 1.0e-5_rk, 'C2 K2 r_north')
        call assert_abs_close(rs, 0.793701_rk, 1.0e-5_rk, 'C2 K2 r_south')
        call poles_s(gb, rn, rs)
        call assert_abs_close(rn, 0.793701_rk, 1.0e-5_rk, 'C2 S2 r_north')
        call profile_s(ga, za, ra, da)
        call profile_s(gb, zb, rb, db)
        call assert_true(maxval(abs(za - zb)) <= 1.0e-15_rk .and. maxval(abs(ra - rb)) <= 1.0e-15_rk &
                .and. maxval(abs(da - db)) <= 1.0e-15_rk, 'C2 the two spheroid routes give one profile')
    end subroutine test_c2_spheroids

    subroutine join_check(params, label)
        real(kind = rk), intent(in) :: params(5)
        character(len = *), intent(in) :: label
        type(geometry_t) :: g
        real(kind = rk) :: r_a, d_a, r_b, d_b
        call resolve_ok(params, g, label)
        call segment_rho_s(g%a_sq(1), 0.0_rk, g%k(1), g%l(1), g%z_join(1), r_a, d_a)
        call segment_rho_s(g%a_sq(3), 0.0_rk, g%k(3), g%l(3), g%z_join(1), r_b, d_b)
        call assert_abs_close(r_a, r_b, 1.0e-13_rk, label // ': rho continuous at z1')
        call assert_abs_close(d_a, d_b, 1.0e-13_rk, label // ': drho/dz continuous at z1')
        call segment_rho_s(g%a_sq(3), 0.0_rk, g%k(3), g%l(3), g%z_join(2), r_a, d_a)
        call segment_rho_s(g%a_sq(2), 0.0_rk, g%k(2), g%l(2), g%z_join(2), r_b, d_b)
        call assert_abs_close(r_a, r_b, 1.0e-13_rk, label // ': rho continuous at z2')
        call assert_abs_close(d_a, d_b, 1.0e-13_rk, label // ': drho/dz continuous at z2')
    end subroutine join_check

    subroutine test_c3_join_continuity()
        call join_check(ASYM_NECKED, 'C3 ASYM_NECKED')
        call join_check(PEAR, 'C3 PEAR')
        call join_check(SPHEROID_MID, 'C3 SPHEROID_MID')
        call join_check(CYLINDER, 'C3 CYLINDER')
    end subroutine test_c3_join_continuity

    subroutine volume_check(params, label)
        real(kind = rk), intent(in) :: params(5)
        character(len = *), intent(in) :: label
        type(geometry_t) :: g
        real(kind = rk) :: vol, zm, q2
        call resolve_ok(params, g, label)
        call gl_body_moments(g, vol, zm, q2)
        call assert_abs_close(vol, 4.0_rk * PI_C / 3.0_rk, 1.0e-11_rk, label // ': V = 4 pi / 3')
        call assert_abs_close(zm, 0.0_rk, 1.0e-11_rk, label // ': z_cm = 0')
    end subroutine volume_check

    subroutine test_c4_volume_and_com()
        call volume_check(SPHERE, 'C4 SPHERE')
        call volume_check(SPHERE_B, 'C4 SPHERE_B')
        call volume_check(SPHEROID_K2, 'C4 SPHEROID_K2')
        call volume_check(SPHEROID_S2, 'C4 SPHEROID_S2')
        call volume_check(SYM_NECKED, 'C4 SYM_NECKED')
        call volume_check(INITIAL, 'C4 INITIAL')
        call volume_check(ASYM_NECKED, 'C4 ASYM_NECKED')
        call volume_check(SPHEROID_MID, 'C4 SPHEROID_MID')
        call volume_check(CYLINDER, 'C4 CYLINDER')
        call volume_check(PEAR, 'C4 PEAR')
        call volume_check(PINCHED, 'C4 PINCHED')
        call volume_check(STAR_FAIL, 'C4 STAR_FAIL')
        call volume_check(DEGENERATE_ONE, 'C4 DEGENERATE_ONE')
        call volume_check(SEPARATED, 'C4 SEPARATED')
        call volume_check(SEPARATED_ASYM, 'C4 SEPARATED_ASYM')
    end subroutine test_c4_volume_and_com

    subroutine derivative_check(params, label)
        real(kind = rk), intent(in) :: params(5)
        character(len = *), intent(in) :: label
        integer(kind = ik), parameter :: N = 1001_ik
        real(kind = rk), parameter :: H = 1.0e-6_rk, KEEP_OUT = 2.0e-5_rk
        type(geometry_t) :: g
        real(kind = rk) :: z(N), rho(N), drho(N), rp, rm, dummy, fd
        integer(kind = ik) :: j, n_checked, n_bad
        call resolve_ok(params, g, label)
        call profile_s(g, z, rho, drho)
        n_checked = 0_ik
        n_bad = 0_ik
        do j = 2_ik, N - 1_ik
            if (min(abs(z(j) - g%z_join(1)), abs(z(j) - g%z_join(2)), abs(z(j) - g%z_south), &
                    abs(z(j) - g%z_north)) <= KEEP_OUT) cycle
            call rho_at_z_s(g, z(j) + H, rp, dummy)
            call rho_at_z_s(g, z(j) - H, rm, dummy)
            fd = (rp - rm) / (2.0_rk * H)
            n_checked = n_checked + 1_ik
            if (abs(drho(j) - fd) > 1.0e-6_rk * (1.0_rk + abs(drho(j)))) n_bad = n_bad + 1_ik
        end do
        call assert_true(n_checked > 900_ik, label // ': most interior nodes checked')
        call assert_true(n_bad == 0_ik, label // ': drho/dz matches central differences')
    end subroutine derivative_check

    subroutine test_c5_profile_derivative()
        call derivative_check(ASYM_NECKED, 'C5 ASYM_NECKED')
        call derivative_check(PEAR, 'C5 PEAR')
    end subroutine test_c5_profile_derivative

    subroutine test_c6_conventions()
        integer(kind = ik), parameter :: N = 1001_ik
        type(geometry_t) :: g
        real(kind = rk) :: z(N), rho(N), drho(N)
        integer(kind = ik) :: j
        logical :: void_zero, outside_positive, any_in_void

        call resolve_ok(SEPARATED_ASYM, g, 'C6 SEPARATED_ASYM')
        call profile_s(g, z, rho, drho)
        void_zero = .true.
        outside_positive = .true.
        any_in_void = .false.
        do j = 2_ik, N - 1_ik
            if (z(j) > g%l(3) - g%w .and. z(j) < g%l(3) + g%w) then
                any_in_void = .true.
                if (.not. (is_zero_f(rho(j)) .and. is_zero_f(drho(j)))) void_zero = .false.
                ! The nodes just outside the void.
                if (.not. (z(j - 1_ik) > g%l(3) - g%w)) outside_positive = outside_positive .and. rho(j - 1_ik) > 0.0_rk
                if (.not. (z(j + 1_ik) < g%l(3) + g%w)) outside_positive = outside_positive .and. rho(j + 1_ik) > 0.0_rk
            end if
        end do
        call assert_true(any_in_void, 'C6 SEPARATED_ASYM: grid nodes inside the void')
        call assert_true(void_zero, 'C6 SEPARATED_ASYM: (0, 0) at every void node')
        call assert_true(outside_positive, 'C6 SEPARATED_ASYM: rho > 0 just outside the void')

        call resolve_ok(PINCHED, g, 'C6 PINCHED')
        call profile_s(g, z, rho, drho)
        call assert_true(is_zero_f(z(501)), 'C6 PINCHED: node 501 is z = 0')
        call assert_true(is_zero_f(rho(501)) .and. is_zero_f(drho(501)), 'C6 PINCHED: (0, 0) at the pinch')
        call assert_true(rho(500) > 0.0_rk .and. rho(502) > 0.0_rk, 'C6 PINCHED: rho > 0 beside the pinch')
    end subroutine test_c6_conventions

end program tqs_param_core_test
