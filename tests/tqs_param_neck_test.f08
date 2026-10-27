!> Neck: the analytic waist of a necked middle, none for convex or flaring
!! middles, rho = 0 at the split plane of pinched and separated shapes.
program tqs_param_neck_test
    use precision_utilities_mod, only: rk
    use tqs_workers_mod, only: geometry_t, neck_s
    use tqs_test_fixtures_mod, only: SYM_NECKED, INITIAL, ASYM_NECKED, SPHEROID_MID, &
            CYLINDER, PEAR, SPHERE, SPHERE_B, SPHEROID_S2, PINCHED, SEPARATED, &
            SEPARATED_ASYM, resolve_ok
    use test_utils_mod, only: assert_true, assert_abs_close, bits_eq_f, is_zero_f, &
            test_summary
    implicit none

    call test_n1_symmetric()
    call test_n2_asymmetric()
    call test_n3_not_found()
    call test_n4_split_plane()

    call test_summary()

contains

    subroutine test_n1_symmetric()
        type(geometry_t) :: g
        real(kind = rk) :: zn, rn
        logical :: found
        call resolve_ok(SYM_NECKED, g, 'N1 SYM_NECKED')
        call neck_s(g, zn, rn, found)
        call assert_true(found, 'N1 SYM_NECKED found')
        call assert_abs_close(zn, 0.0_rk, 1.0e-15_rk, 'N1 SYM_NECKED z_neck = 0')
        call assert_true(bits_eq_f(rn, sqrt(g%a_sq(3))), 'N1 SYM_NECKED rho_neck = a3')
        call resolve_ok(INITIAL, g, 'N1 INITIAL')
        call neck_s(g, zn, rn, found)
        call assert_true(found, 'N1 INITIAL found')
        call assert_abs_close(zn, 0.0_rk, 1.0e-15_rk, 'N1 INITIAL z_neck = 0')
        call assert_true(bits_eq_f(rn, sqrt(g%a_sq(3))), 'N1 INITIAL rho_neck = a3')
    end subroutine test_n1_symmetric

    subroutine test_n2_asymmetric()
        type(geometry_t) :: g
        real(kind = rk) :: zn, rn
        logical :: found
        call resolve_ok(ASYM_NECKED, g, 'N2 ASYM_NECKED')
        call neck_s(g, zn, rn, found)
        call assert_true(found, 'N2 found')
        call assert_true(bits_eq_f(zn, g%l(3)), 'N2 z_neck = l3')
        call assert_abs_close(rn, 0.674685_rk, 1.0e-5_rk, 'N2 rho_neck')
    end subroutine test_n2_asymmetric

    subroutine expect_none(params, label)
        real(kind = rk), intent(in) :: params(5)
        character(len = *), intent(in) :: label
        type(geometry_t) :: g
        real(kind = rk) :: zn, rn
        logical :: found
        call resolve_ok(params, g, label)
        call neck_s(g, zn, rn, found)
        call assert_true(.not. found .and. is_zero_f(zn) .and. is_zero_f(rn), &
                label // ': not found, zero outputs')
    end subroutine expect_none

    subroutine test_n3_not_found()
        call expect_none(SPHEROID_MID, 'N3 SPHEROID_MID')
        call expect_none(CYLINDER, 'N3 CYLINDER')
        call expect_none(PEAR, 'N3 PEAR')
        call expect_none(SPHERE, 'N3 SPHERE')
        call expect_none(SPHERE_B, 'N3 SPHERE_B')
        call expect_none(SPHEROID_S2, 'N3 SPHEROID_S2')
    end subroutine test_n3_not_found

    subroutine test_n4_split_plane()
        type(geometry_t) :: g
        real(kind = rk) :: zn, rn
        logical :: found
        call resolve_ok(PINCHED, g, 'N4 PINCHED')
        call neck_s(g, zn, rn, found)
        call assert_true(found .and. is_zero_f(rn), 'N4 PINCHED found with rho_neck = 0')
        call assert_abs_close(zn, 0.0_rk, 1.0e-15_rk, 'N4 PINCHED z_neck = 0')
        call resolve_ok(SEPARATED, g, 'N4 SEPARATED')
        call neck_s(g, zn, rn, found)
        call assert_true(found .and. is_zero_f(rn), 'N4 SEPARATED found with rho_neck = 0')
        call assert_abs_close(zn, 0.0_rk, 1.0e-15_rk, 'N4 SEPARATED z_neck = 0')
        call resolve_ok(SEPARATED_ASYM, g, 'N4 SEPARATED_ASYM')
        call neck_s(g, zn, rn, found)
        call assert_true(found .and. is_zero_f(rn), 'N4 SEPARATED_ASYM found with rho_neck = 0')
        call assert_true(bits_eq_f(zn, g%l(3)), 'N4 SEPARATED_ASYM z_neck = l3')
    end subroutine test_n4_split_plane

end program tqs_param_neck_test
