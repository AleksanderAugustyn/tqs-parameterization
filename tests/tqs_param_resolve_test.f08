!> Resolve: coordinates -> geometry. Part 1 is the unit stage (u = 1,
!! unshifted): validation, joins, the root rule, the degeneracies and the
!! classification.
program tqs_param_resolve_test
    use precision_utilities_mod, only: ik, rk
    use tqs_workers_mod, only: geometry_t, TQS_MID_SPHEROID, TQS_MID_CYLINDER, &
            TQS_MID_ONE_SHEET, TQS_MID_CONE, TQS_MID_TWO_SHEET, &
            TQS_ERROR_INVALID_COORDINATES, TQS_ERROR_NO_GEOMETRY
    use tqs_resolve_mod, only: resolve_unit_s
    use test_utils_mod, only: assert_true, assert_int_eq, assert_abs_close, &
            arrays_bits_eq_f, bits_eq_f, is_zero_f, test_summary
    implicit none

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
    real(kind = rk), parameter :: STAR_FAIL(5) = [3.0_rk, -0.6_rk, 1.0_rk, 0.9_rk, 0.0_rk]
    real(kind = rk), parameter :: DEGENERATE_ONE(5) = [0.5_rk, 2.0_rk, 1.5_rk, -0.5_rk, 1.0_rk]
    real(kind = rk), parameter :: BAD_K(5) = [0.0_rk, 1.0_rk, 0.5_rk, 0.0_rk, 1.0_rk]
    real(kind = rk), parameter :: BAD_A2(5) = [0.5_rk, -0.6_rk, 1.0_rk, 2.0_rk, 0.0_rk]
    real(kind = rk), parameter :: BAD_S1(5) = [-0.1_rk, -0.6_rk, 1.0_rk, 0.0_rk, 0.0_rk]
    real(kind = rk), parameter :: NO_GEOM(5) = [0.3_rk, 1.5_rk, 1.5_rk, 0.0_rk, 0.0_rk]
    real(kind = rk), parameter :: NO_CYL(5) = [0.3_rk, 0.0_rk, 1.0_rk, 0.2_rk, 0.0_rk]

    call test_r1_symmetric()
    call test_r2_pear()
    call test_r3_degenerate_one()
    call test_r4_both_degenerate()
    call test_r5_cylinder()
    call test_r6_sigma1_zero()
    call test_r7_invalid_coordinates()
    call test_r8_classification()
    call test_r9_determinism()

    call test_summary()

contains

    !> Resolve at the unit stage and assert the status.
    subroutine unit_ok(params, geom, label)
        real(kind = rk), intent(in) :: params(5)
        type(geometry_t), intent(out) :: geom
        character(len = *), intent(in) :: label
        integer(kind = ik) :: status
        call resolve_unit_s(params, geom, status)
        call assert_int_eq(status, 0_ik, label // ': status 0')
    end subroutine unit_ok

    subroutine expect_status(params, want, label)
        real(kind = rk), intent(in) :: params(5)
        integer(kind = ik), intent(in) :: want
        character(len = *), intent(in) :: label
        type(geometry_t) :: geom
        integer(kind = ik) :: status
        call resolve_unit_s(params, geom, status)
        call assert_int_eq(status, want, label)
    end subroutine expect_status

    subroutine test_r1_symmetric()
        type(geometry_t) :: geom
        call unit_ok(SYM_NECKED, geom, 'R1 SYM_NECKED')
        call assert_true(is_zero_f(geom%l(3)), 'R1 SYM_NECKED l3 = 0 exactly')
        call assert_abs_close(geom%z_join(1), -geom%z_join(2), 1.0e-15_rk, 'R1 SYM_NECKED z1 = -z2')
        call unit_ok(INITIAL, geom, 'R1 INITIAL')
        call assert_true(is_zero_f(geom%l(3)), 'R1 INITIAL l3 = 0 exactly')
        call assert_abs_close(geom%z_join(1), -geom%z_join(2), 1.0e-15_rk, 'R1 INITIAL z1 = -z2')
        call unit_ok(SPHERE_B, geom, 'R1 SPHERE_B')
        call assert_true(is_zero_f(geom%l(3)), 'R1 SPHERE_B l3 = 0 exactly')
        call assert_abs_close(geom%z_join(1), -geom%z_join(2), 1.0e-15_rk, 'R1 SPHERE_B z1 = -z2')
    end subroutine test_r1_symmetric

    subroutine test_r2_pear()
        type(geometry_t) :: geom
        call unit_ok(PEAR, geom, 'R2 PEAR')
        call assert_abs_close(geom%l(3), -1.4_rk, 1.0e-12_rk, 'R2 PEAR l3')
        call assert_abs_close(geom%z_join(1), -0.95_rk, 1.0e-12_rk, 'R2 PEAR z1')
        call assert_abs_close(geom%z_join(2), -0.45_rk, 1.0e-12_rk, 'R2 PEAR z2')
        call assert_abs_close(geom%a_sq(3), -0.105_rk, 1.0e-12_rk, 'R2 PEAR a3^2')
        call assert_int_eq(geom%mid_kind, TQS_MID_TWO_SHEET, 'R2 PEAR two_sheet')
        call assert_true(.not. geom%waist_inside, 'R2 PEAR waist outside')
        call assert_true(.not. geom%is_separated, 'R2 PEAR not separated')
        call assert_abs_close(geom%u, 1.0_rk, 0.0_rk, 'R2 PEAR unit stage u = 1')
    end subroutine test_r2_pear

    subroutine test_r3_degenerate_one()
        type(geometry_t) :: geom
        integer(kind = ik) :: status, n_valid
        call resolve_unit_s(DEGENERATE_ONE, geom, status, n_valid)
        call assert_int_eq(status, 0_ik, 'R3 DEGENERATE_ONE status 0')
        call assert_abs_close(geom%z_join(1), -0.75_rk, 1.0e-12_rk, 'R3 z1 = -0.75')
        call assert_true(bits_eq_f(geom%z_join(1), geom%z_join(2)), 'R3 z1 = z2 bitwise')
        call assert_abs_close(geom%l(3), -0.25_rk, 1.0e-12_rk, 'R3 l3 = -0.25')
        call assert_abs_close(geom%a_sq(3), 0.75_rk, 1.0e-12_rk, 'R3 a3^2 = 0.75')
        call assert_int_eq(geom%mid_kind, TQS_MID_SPHEROID, 'R3 spheroid middle')
        call assert_int_eq(n_valid, 0_ik, 'R3 degenerate path: n_valid_roots = 0')
    end subroutine test_r3_degenerate_one

    subroutine test_r4_both_degenerate()
        type(geometry_t) :: geom
        call unit_ok(SPHERE, geom, 'R4 SPHERE')
        call assert_abs_close(geom%z_join(1), 0.0_rk, 1.0e-15_rk, 'R4 SPHERE z1 = 0')
        call assert_abs_close(geom%z_join(2), 0.0_rk, 1.0e-15_rk, 'R4 SPHERE z2 = 0')
        call assert_abs_close(geom%l(3), 0.0_rk, 1.0e-15_rk, 'R4 SPHERE l3 = 0')
        call assert_abs_close(geom%a_sq(3), 1.0_rk, 1.0e-15_rk, 'R4 SPHERE a3^2 = 1')
        call assert_abs_close(geom%k(3), 1.0_rk, 0.0_rk, 'R4 SPHERE k3 = sigma2')
        call unit_ok(SPHEROID_K2, geom, 'R4 SPHEROID_K2')
        call assert_abs_close(geom%z_join(1), 0.0_rk, 1.0e-15_rk, 'R4 SPHEROID_K2 z1 = 0')
        call assert_abs_close(geom%z_join(2), 0.0_rk, 1.0e-15_rk, 'R4 SPHEROID_K2 z2 = 0')
        call assert_abs_close(geom%l(3), 0.0_rk, 1.0e-15_rk, 'R4 SPHEROID_K2 l3 = 0')
        call assert_abs_close(geom%a_sq(3), 1.0_rk, 1.0e-15_rk, 'R4 SPHEROID_K2 a3^2 = 1')
        call expect_status(NO_GEOM, TQS_ERROR_NO_GEOMETRY, 'R4 NO_GEOM -> 101')
        call expect_status([0.0_rk, 1.5_rk, 1.5_rk, 0.2_rk, 0.0_rk], TQS_ERROR_NO_GEOMETRY, &
                'R4 both degenerate with alpha2 /= 0 -> 101')
    end subroutine test_r4_both_degenerate

    subroutine test_r5_cylinder()
        type(geometry_t) :: geom
        call unit_ok(CYLINDER, geom, 'R5 CYLINDER')
        call assert_int_eq(geom%mid_kind, TQS_MID_CYLINDER, 'R5 cylinder kind')
        call assert_abs_close(geom%z_join(1), -0.3_rk, 1.0e-15_rk, 'R5 z1 = -0.3')
        call assert_abs_close(geom%z_join(2), 0.3_rk, 1.0e-15_rk, 'R5 z2 = 0.3')
        call assert_abs_close(geom%a_sq(3), 1.0_rk, 1.0e-15_rk, 'R5 a3^2 = 1')
        call expect_status(NO_CYL, TQS_ERROR_NO_GEOMETRY, 'R5 NO_CYL -> 101')
    end subroutine test_r5_cylinder

    subroutine test_r6_sigma1_zero()
        type(geometry_t) :: geom
        call unit_ok(SPHEROID_S2, geom, 'R6 SPHEROID_S2')
        call assert_abs_close(geom%z_join(1), 0.0_rk, 1.0e-15_rk, 'R6 z1 = 0')
        call assert_abs_close(geom%z_join(2), 0.0_rk, 1.0e-15_rk, 'R6 z2 = 0')
        call assert_abs_close(geom%l(3), 0.0_rk, 1.0e-15_rk, 'R6 l3 = 0')
        call assert_abs_close(geom%a_sq(3), 1.0_rk, 1.0e-15_rk, 'R6 a3^2 = 1')
        call assert_int_eq(geom%mid_kind, TQS_MID_ONE_SHEET, 'R6 one_sheet')
        call assert_true(.not. geom%waist_inside, 'R6 waist not inside a zero-length middle')
    end subroutine test_r6_sigma1_zero

    subroutine test_r7_invalid_coordinates()
        call expect_status(BAD_K, TQS_ERROR_INVALID_COORDINATES, 'R7 BAD_K -> 100')
        call expect_status(BAD_A2, TQS_ERROR_INVALID_COORDINATES, 'R7 BAD_A2 -> 100')
        call expect_status(BAD_S1, TQS_ERROR_INVALID_COORDINATES, 'R7 BAD_S1 -> 100')
        call expect_status([0.5_rk, -0.6_rk, 1.0_rk, -2.0_rk, 0.0_rk], &
                TQS_ERROR_INVALID_COORDINATES, 'R7 alpha2 = -2 -> 100')
    end subroutine test_r7_invalid_coordinates

    subroutine test_r8_classification()
        type(geometry_t) :: geom
        call unit_ok(SPHEROID_MID, geom, 'R8 SPHEROID_MID')
        call assert_int_eq(geom%mid_kind, TQS_MID_SPHEROID, 'R8 SPHEROID_MID spheroid')
        call unit_ok(SYM_NECKED, geom, 'R8 SYM_NECKED')
        call assert_int_eq(geom%mid_kind, TQS_MID_ONE_SHEET, 'R8 SYM_NECKED one_sheet')
        call assert_true(geom%waist_inside, 'R8 SYM_NECKED waist inside')
        call unit_ok(PINCHED, geom, 'R8 PINCHED')
        call assert_int_eq(geom%mid_kind, TQS_MID_CONE, 'R8 PINCHED cone')
        call assert_true(geom%is_pinched, 'R8 PINCHED is_pinched')
        call assert_true(.not. geom%is_separated, 'R8 PINCHED not separated')
        call unit_ok(SEPARATED, geom, 'R8 SEPARATED')
        call assert_int_eq(geom%mid_kind, TQS_MID_TWO_SHEET, 'R8 SEPARATED two_sheet')
        call assert_true(geom%is_separated, 'R8 SEPARATED is_separated')
        call assert_true(geom%z_join(1) < geom%l(3) - geom%w .and. &
                geom%l(3) + geom%w < geom%z_join(2), 'R8 SEPARATED void strictly inside the middle')
        call unit_ok(STAR_FAIL, geom, 'R8 STAR_FAIL')
        call assert_int_eq(geom%mid_kind, TQS_MID_ONE_SHEET, 'R8 STAR_FAIL one_sheet')
        call assert_true(geom%waist_inside, 'R8 STAR_FAIL waist inside')
        call assert_true(.not. geom%is_separated, 'R8 STAR_FAIL not separated')
    end subroutine test_r8_classification

    subroutine test_r9_determinism()
        type(geometry_t) :: g1, g2
        call unit_ok(ASYM_NECKED, g1, 'R9 first')
        call unit_ok(ASYM_NECKED, g2, 'R9 second')
        call assert_true(geometries_bits_eq(g1, g2), 'R9 two calls bitwise equal')
    end subroutine test_r9_determinism

    logical function geometries_bits_eq(g1, g2) result(same)
        type(geometry_t), intent(in) :: g1, g2
        same = bits_eq_f(g1%u, g2%u) .and. arrays_bits_eq_f(g1%a_sq, g2%a_sq) &
                .and. arrays_bits_eq_f(g1%k, g2%k) .and. arrays_bits_eq_f(g1%c, g2%c) &
                .and. arrays_bits_eq_f(g1%l, g2%l) .and. arrays_bits_eq_f(g1%z_join, g2%z_join) &
                .and. bits_eq_f(g1%z_south, g2%z_south) .and. bits_eq_f(g1%z_north, g2%z_north) &
                .and. bits_eq_f(g1%w, g2%w) .and. g1%mid_kind == g2%mid_kind &
                .and. (g1%waist_inside .eqv. g2%waist_inside) &
                .and. (g1%is_separated .eqv. g2%is_separated) &
                .and. (g1%is_pinched .eqv. g2%is_pinched)
    end function geometries_bits_eq

end program tqs_param_resolve_test
