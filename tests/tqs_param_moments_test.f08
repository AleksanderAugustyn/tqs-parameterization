!> Moller coordinates: q2 against quadrature, the eps round trip, alpha_g
!! against the completed end-spheroid volumes, and the d conventions.
program tqs_param_moments_test
    use precision_utilities_mod, only: rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use tqs_workers_mod, only: geometry_t, moller_coordinates_s, rho_at_z_s
    use tqs_test_fixtures_mod, only: SPHERE, SPHERE_B, SPHEROID_K2, SPHEROID_S2, &
            SYM_NECKED, INITIAL, ASYM_NECKED, SPHEROID_MID, CYLINDER, PEAR, PINCHED, &
            SEPARATED, SEPARATED_ASYM, STAR_FAIL, DEGENERATE_ONE, pv, resolve_ok, &
            gl_body_moments
    use test_utils_mod, only: assert_true, assert_abs_close, bits_eq_f, is_zero_f, &
            test_summary
    implicit none

    call test_m1_q2()
    call test_m2_eps_round_trip()
    call test_m3_alpha_g()
    call test_m4_d()

    call test_summary()

contains

    subroutine q2_check(params, label, q2_out)
        real(kind = rk), intent(in) :: params(5)
        character(len = *), intent(in) :: label
        real(kind = rk), intent(out), optional :: q2_out
        type(geometry_t) :: g
        real(kind = rk) :: q2, d, ag, e1, e2, vol, zm, q2_gl
        call resolve_ok(params, g, label)
        call moller_coordinates_s(g, q2, d, ag, e1, e2)
        call gl_body_moments(g, vol, zm, q2_gl)
        call assert_abs_close(q2, q2_gl, 1.0e-12_rk, label // ': q2 = quadrature')
        if (present(q2_out)) q2_out = q2
    end subroutine q2_check

    subroutine test_m1_q2()
        real(kind = rk) :: q2
        call q2_check(SPHERE, 'M1 SPHERE')
        call q2_check(SPHERE_B, 'M1 SPHERE_B')
        call q2_check(SPHEROID_K2, 'M1 SPHEROID_K2')
        call q2_check(SPHEROID_S2, 'M1 SPHEROID_S2')
        call q2_check(SYM_NECKED, 'M1 SYM_NECKED', q2)
        call assert_abs_close(q2, 0.380550_rk, 1.0e-5_rk, 'M1 SYM_NECKED cross-check')
        call q2_check(INITIAL, 'M1 INITIAL')
        call q2_check(ASYM_NECKED, 'M1 ASYM_NECKED', q2)
        call assert_abs_close(q2, 0.364920_rk, 1.0e-5_rk, 'M1 ASYM_NECKED cross-check')
        call q2_check(SPHEROID_MID, 'M1 SPHEROID_MID')
        call q2_check(CYLINDER, 'M1 CYLINDER')
        call q2_check(PEAR, 'M1 PEAR')
        call q2_check(PINCHED, 'M1 PINCHED')
        call q2_check(SEPARATED, 'M1 SEPARATED')
        call q2_check(SEPARATED_ASYM, 'M1 SEPARATED_ASYM')
        call q2_check(STAR_FAIL, 'M1 STAR_FAIL')
        call q2_check(DEGENERATE_ONE, 'M1 DEGENERATE_ONE')
    end subroutine test_m1_q2

    subroutine test_m2_eps_round_trip()
        real(kind = rk), parameter :: EPS1 = 0.3_rk, EPS2 = -0.1_rk
        real(kind = rk) :: r1, r2, k1, k2, q2, d, ag, e1, e2, p(5)
        type(geometry_t) :: g
        r1 = (3.0_rk - 2.0_rk * EPS1) / (3.0_rk + EPS1)
        r2 = (3.0_rk - 2.0_rk * EPS2) / (3.0_rk + EPS2)
        k1 = r1 * r1
        k2 = r2 * r2
        p = pv(1.6_rk, -0.6_rk, 0.5_rk * (k1 + k2), 0.0_rk, k1 - k2)
        call resolve_ok(p, g, 'M2 eps fixture')
        call moller_coordinates_s(g, q2, d, ag, e1, e2)
        call assert_abs_close(e1, EPS1, 1.0e-13_rk, 'M2 eps_f1 recovered')
        call assert_abs_close(e2, EPS2, 1.0e-13_rk, 'M2 eps_f2 recovered')
    end subroutine test_m2_eps_round_trip

    subroutine alpha_g_check(params, label, symmetric)
        real(kind = rk), intent(in) :: params(5)
        character(len = *), intent(in) :: label
        logical, intent(in) :: symmetric
        type(geometry_t) :: g
        real(kind = rk) :: q2, d, ag, e1, e2, v1, v2
        call resolve_ok(params, g, label)
        call moller_coordinates_s(g, q2, d, ag, e1, e2)
        v1 = 4.0_rk * PI_C / 3.0_rk * g%a_sq(1) * g%c(1)
        v2 = 4.0_rk * PI_C / 3.0_rk * g%a_sq(2) * g%c(2)
        if (symmetric) then
            call assert_abs_close(ag, 0.0_rk, 1.0e-15_rk, label // ': alpha_g = 0')
        else
            call assert_abs_close(ag, (v1 - v2) / (v1 + v2), 1.0e-14_rk, label // ': alpha_g from end volumes')
        end if
    end subroutine alpha_g_check

    subroutine test_m3_alpha_g()
        call alpha_g_check(ASYM_NECKED, 'M3 ASYM_NECKED', .false.)
        call alpha_g_check(PEAR, 'M3 PEAR', .false.)
        call alpha_g_check(SYM_NECKED, 'M3 SYM_NECKED', .true.)
        call alpha_g_check(INITIAL, 'M3 INITIAL', .true.)
        call alpha_g_check(SPHERE, 'M3 SPHERE', .true.)
    end subroutine test_m3_alpha_g

    subroutine test_m4_d()
        type(geometry_t) :: g
        real(kind = rk) :: q2, d, ag, e1, e2, r1, r2, dummy

        call resolve_ok(SYM_NECKED, g, 'M4 SYM_NECKED')
        call moller_coordinates_s(g, q2, d, ag, e1, e2)
        call assert_true(bits_eq_f(d, 2.0_rk * sqrt(g%a_sq(3))), 'M4 SYM_NECKED d = 2 a3')

        call resolve_ok(PINCHED, g, 'M4 PINCHED')
        call moller_coordinates_s(g, q2, d, ag, e1, e2)
        call assert_true(is_zero_f(d), 'M4 PINCHED d = 0')
        call resolve_ok(SEPARATED, g, 'M4 SEPARATED')
        call moller_coordinates_s(g, q2, d, ag, e1, e2)
        call assert_true(is_zero_f(d), 'M4 SEPARATED d = 0')

        call resolve_ok(SPHEROID_MID, g, 'M4 SPHEROID_MID')
        call moller_coordinates_s(g, q2, d, ag, e1, e2)
        call rho_at_z_s(g, g%z_join(1), r1, dummy)
        call rho_at_z_s(g, g%z_join(2), r2, dummy)
        call assert_true(bits_eq_f(d, 2.0_rk * min(r1, r2)), 'M4 SPHEROID_MID d = narrowest join')
        call resolve_ok(PEAR, g, 'M4 PEAR')
        call moller_coordinates_s(g, q2, d, ag, e1, e2)
        call rho_at_z_s(g, g%z_join(1), r1, dummy)
        call rho_at_z_s(g, g%z_join(2), r2, dummy)
        call assert_true(bits_eq_f(d, 2.0_rk * min(r1, r2)), 'M4 PEAR d = narrowest join')

        call resolve_ok(SPHERE, g, 'M4 SPHERE')
        call moller_coordinates_s(g, q2, d, ag, e1, e2)
        call assert_abs_close(d, 2.0_rk, 1.0e-14_rk, 'M4 SPHERE d = 2')
    end subroutine test_m4_d

end program tqs_param_moments_test
