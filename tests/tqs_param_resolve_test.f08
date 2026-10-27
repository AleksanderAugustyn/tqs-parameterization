!> Resolve: coordinates -> geometry. Part 1 is the unit stage (u = 1,
!! unshifted): validation, joins, the root rule, the degeneracies and the
!! classification. Part 2 is the R0 stage (volume scale, centre of mass,
!! poles) and the scission map.
program tqs_param_resolve_test
    use precision_utilities_mod, only: ik, rk
    use tqs_workers_mod, only: geometry_t, TQS_MID_SPHEROID, TQS_MID_CYLINDER, &
            TQS_MID_ONE_SHEET, TQS_MID_CONE, TQS_MID_TWO_SHEET, &
            TQS_ERROR_INVALID_COORDINATES, TQS_ERROR_NO_GEOMETRY, TQS_PINCH_TOL
    use tqs_resolve_mod, only: resolve_unit_s, resolve_s, alpha1_f, scission_sigma1_s
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
    call test_r10_scission_anchors()
    call test_r11_locus()
    call test_r12_star_above_t()
    call test_r13_eq16()
    call test_r14_continuity_walks()
    call test_r15_r0_cross_checks()
    call test_r16_symmetric_after_shift()

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
        ! waist_inside is a sigma2 < 0 notion: a convex middle has no waist.
        call assert_true(.not. geom%waist_inside, 'R8 SPHEROID_MID waist_inside false')
        call unit_ok(CYLINDER, geom, 'R8 CYLINDER')
        call assert_true(.not. geom%waist_inside, 'R8 CYLINDER waist_inside false')
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

    !---------------------------------------------------------------------------
    ! Part 2: R0 stage and the scission map
    !---------------------------------------------------------------------------

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

    subroutine resolve_ok(params, geom, label)
        real(kind = rk), intent(in) :: params(5)
        type(geometry_t), intent(out) :: geom
        character(len = *), intent(in) :: label
        integer(kind = ik) :: status
        call resolve_s(params, geom, status)
        call assert_int_eq(status, 0_ik, label // ': status 0')
    end subroutine resolve_ok

    subroutine test_r10_scission_anchors()
        real(kind = rk) :: t, star, t7, star7
        logical :: defined, defined7
        integer(kind = ik) :: status

        call scission_sigma1_s([1.0_rk, -1.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], t, star, defined, status)
        call assert_int_eq(status, 0_ik, 'R10 sphere ends status 0')
        call assert_abs_close(t, 2.0_rk, 1.0e-14_rk, 'R10 sphere ends sigma1_t = 2')
        call assert_abs_close(star, 2.0_rk * sqrt(2.0_rk), 1.0e-14_rk, 'R10 sphere ends sigma1_star = 2 sqrt 2')
        call assert_true(defined, 'R10 sphere ends defined')

        call scission_sigma1_s([1.0_rk, -0.6_rk, 1.0_rk, 0.3_rk, 0.2_rk], t, star, defined, status)
        call assert_abs_close(t, 1.994300_rk, 1.0e-5_rk, 'R10 ASYM ends sigma1_t')
        call assert_abs_close(star, 3.257672_rk, 1.0e-5_rk, 'R10 ASYM ends sigma1_star')

        call scission_sigma1_s([1.0_rk, 0.5_rk, 1.0_rk, 0.3_rk, 0.2_rk], t, star, defined, status)
        call assert_int_eq(status, 0_ik, 'R10 sigma2 > 0 status 0')
        call assert_true(.not. defined, 'R10 sigma2 > 0 not defined')
        call assert_true(is_zero_f(star), 'R10 sigma2 > 0 sigma1_star = 0')

        call scission_sigma1_s(BAD_K, t, star, defined, status)
        call assert_int_eq(status, TQS_ERROR_INVALID_COORDINATES, 'R10 BAD_K -> 100')
        call assert_true(is_zero_f(t) .and. is_zero_f(star) .and. .not. defined, &
                'R10 BAD_K zero-filled outputs')

        call scission_sigma1_s([1.0_rk, -0.6_rk, 1.0_rk, 0.3_rk, 0.2_rk], t, star, defined, status)
        call scission_sigma1_s([7.0_rk, -0.6_rk, 1.0_rk, 0.3_rk, 0.2_rk], t7, star7, defined7, status)
        call assert_true(bits_eq_f(t, t7) .and. bits_eq_f(star, star7) .and. (defined .eqv. defined7), &
                'R10 sigma1 does not enter the scission map')
    end subroutine test_r10_scission_anchors

    subroutine test_r11_locus()
        real(kind = rk) :: p(5)
        real(kind = rk) :: t, star
        logical :: defined
        integer(kind = ik) :: status
        type(geometry_t) :: geom

        call scission_sigma1_s([0.0_rk, -0.6_rk, 1.0_rk, 0.3_rk, 0.2_rk], t, star, defined, status)
        p = pv(0.98_rk * star, -0.6_rk, 1.0_rk, 0.3_rk, 0.2_rk)
        call resolve_ok(p, geom, 'R11 below')
        call assert_int_eq(geom%mid_kind, TQS_MID_ONE_SHEET, 'R11 below: one_sheet')
        call assert_true(geom%waist_inside, 'R11 below: waist inside')
        p = pv(star, -0.6_rk, 1.0_rk, 0.3_rk, 0.2_rk)
        call resolve_ok(p, geom, 'R11 at')
        call assert_int_eq(geom%mid_kind, TQS_MID_CONE, 'R11 at: cone')
        call assert_true(geom%is_pinched, 'R11 at: pinched')
        call assert_true(abs(geom%a_sq(3)) <= TQS_PINCH_TOL * geom%u**2, 'R11 at: |a3^2| within the pinch tol')
        p = pv(1.02_rk * star, -0.6_rk, 1.0_rk, 0.3_rk, 0.2_rk)
        call resolve_ok(p, geom, 'R11 above')
        call assert_int_eq(geom%mid_kind, TQS_MID_TWO_SHEET, 'R11 above: two_sheet')
        call assert_true(geom%is_separated, 'R11 above: separated')
    end subroutine test_r11_locus

    subroutine test_r12_star_above_t()
        real(kind = rk) :: p(5)
        real(kind = rk), parameter :: S2_SWEEP(5) = [-0.1_rk, -0.5_rk, -1.0_rk, -5.0_rk, -100.0_rk]
        real(kind = rk) :: t, star
        logical :: defined
        integer(kind = ik) :: status, j
        type(geometry_t) :: geom
        character(len = 32) :: tag

        do j = 1_ik, size(S2_SWEEP, kind = ik)
            write(tag, '(A,ES9.2)') 'R12 s2 = ', S2_SWEEP(j)
            p = pv(0.0_rk, S2_SWEEP(j), 1.0_rk, 0.3_rk, 0.2_rk)
            call scission_sigma1_s(p, t, star, defined, status)
            call assert_true(star > t, trim(tag) // ': sigma1_star > sigma1_t')
        end do

        call scission_sigma1_s([0.0_rk, -1.0e6_rk, 1.0_rk, 0.3_rk, 0.2_rk], t, star, defined, status)
        call assert_abs_close(star, t, 1.0e-5_rk, 'R12 s2 = -1e6: sigma1_star -> sigma1_t')
        p = pv(star, -1.0e6_rk, 1.0_rk, 0.3_rk, 0.2_rk)
        call resolve_ok(p, geom, 'R12 s2 = -1e6 at the locus')
        call assert_abs_close(geom%z_join(1), geom%l(1) + geom%c(1), 1.0e-5_rk, &
                'R12 s2 = -1e6: join 1 at the touching pole')
        call assert_abs_close(geom%z_join(2), geom%l(2) - geom%c(2), 1.0e-5_rk, &
                'R12 s2 = -1e6: join 2 at the touching pole')
    end subroutine test_r12_star_above_t

    subroutine test_r13_eq16()
        real(kind = rk) :: p(5)
        real(kind = rk), parameter :: S3 = 1.0_rk, A2 = 0.3_rk, A3 = 0.2_rk
        real(kind = rk) :: num, eq16, t, star
        logical :: defined
        integer(kind = ik) :: status

        num = (A2 * A3 - 2.0_rk * A2 * S3 + 2.0_rk * A3 - 4.0_rk * S3) &
                / (A2 * A3 + 2.0_rk * A2 * S3 - 2.0_rk * A3 - 4.0_rk * S3)
        eq16 = (sqrt(num) + 1.0_rk) / sqrt((A3 - 2.0_rk * S3) / (A2 - 2.0_rk))
        p = pv(0.0_rk, -0.6_rk, S3, A2, A3)
        call scission_sigma1_s(p, t, star, defined, status)
        call assert_abs_close(eq16, t, 1.0e-12_rk, 'R13 Moller eq. 16 as printed = c1 + c2')
    end subroutine test_r13_eq16

    !> Walk params(slot) over first + step * (0 .. n - 1); every point valid, the
    !! per-step change of z1, z2, l3, a3^2 and u bounded by max_step. Returns the
    !! kinds seen and whether the waist was ever inside.
    subroutine walk(base, slot, first, step, n, max_step, label, kinds_seen, any_inside, kinds)
        real(kind = rk), intent(in) :: base(5), first, step, max_step
        integer(kind = ik), intent(in) :: slot, n
        character(len = *), intent(in) :: label
        logical, intent(out) :: kinds_seen(5), any_inside
        integer(kind = ik), intent(out) :: kinds(n)
        real(kind = rk) :: p(5), prev(5), cur(5), worst
        type(geometry_t) :: geom
        integer(kind = ik) :: j, status, n_bad
        kinds_seen = .false.
        any_inside = .false.
        n_bad = 0_ik
        worst = 0.0_rk
        prev = 0.0_rk
        do j = 1_ik, n
            p = base
            p(slot) = first + step * real(j - 1_ik, rk)
            call resolve_s(p, geom, status)
            if (status /= 0_ik) then
                n_bad = n_bad + 1_ik
                kinds(j) = 0_ik
                cycle
            end if
            kinds(j) = geom%mid_kind
            kinds_seen(geom%mid_kind) = .true.
            any_inside = any_inside .or. geom%waist_inside
            cur(1:2) = geom%z_join
            cur(3) = geom%l(3)
            cur(4) = geom%a_sq(3)
            cur(5) = geom%u
            if (j > 1_ik) worst = max(worst, maxval(abs(cur - prev)))
            prev = cur
        end do
        call assert_int_eq(n_bad, 0_ik, label // ': every point valid')
        call assert_true(worst <= max_step, label // ': per-step change bounded')
    end subroutine walk

    subroutine test_r14_continuity_walks()
        integer(kind = ik), parameter :: N1 = 81_ik, N2 = 61_ik, N3 = 56_ik, N4 = 81_ik
        logical :: seen(5), inside
        integer(kind = ik) :: k1(N1), k2(N2), k3(N3), k4(N4), j, first_two

        ! sigma1 through the scission locus: one_sheet, then two_sheet, never back.
        call walk([0.0_rk, -0.6_rk, 1.0_rk, 0.3_rk, 0.2_rk], 1_ik, 2.8_rk, 0.01_rk, N1, 0.02_rk, &
                'R14 sigma1 walk', seen, inside, k1)
        call assert_int_eq(k1(1), TQS_MID_ONE_SHEET, 'R14 sigma1 walk starts one_sheet')
        call assert_int_eq(k1(N1), TQS_MID_TWO_SHEET, 'R14 sigma1 walk ends two_sheet')
        first_two = max(1_ik, findloc(k1, TQS_MID_TWO_SHEET, dim = 1, kind = ik))
        call assert_true(all(k1(first_two:) == TQS_MID_TWO_SHEET), 'R14 sigma1 walk never returns to one_sheet')
        call assert_true(all(k1(:first_two - 1_ik) == TQS_MID_ONE_SHEET .or. k1(:first_two - 1_ik) == TQS_MID_CONE), &
                'R14 sigma1 walk one_sheet before the locus')

        ! sigma2 through 0 on an integer grid, so sigma2 = 0 is hit exactly.
        call walk([0.6_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.3_rk], 2_ik, -0.3_rk, 0.01_rk, N2, 0.02_rk, &
                'R14 sigma2 walk', seen, inside, k2)
        do j = 1_ik, N2
            if (j == 31_ik) then
                call assert_int_eq(k2(j), TQS_MID_CYLINDER, 'R14 sigma2 walk: cylinder at 0')
            else if (j < 31_ik) then
                call assert_int_eq(k2(j), TQS_MID_ONE_SHEET, 'R14 sigma2 walk: one_sheet below 0')
            else
                call assert_int_eq(k2(j), TQS_MID_SPHEROID, 'R14 sigma2 walk: spheroid above 0')
            end if
        end do

        ! alpha2 along the pear family: a3^2 crosses 0 with the waist outside.
        call walk([1.0_rk, -1.0_rk, 1.0_rk, 0.0_rk, 0.0_rk], 4_ik, -1.45_rk, 0.01_rk, N3, 0.03_rk, &
                'R14 alpha2 walk', seen, inside, k3)
        call assert_true(.not. inside, 'R14 alpha2 walk: waist outside throughout')
        call assert_true(seen(TQS_MID_ONE_SHEET) .and. seen(TQS_MID_TWO_SHEET), &
                'R14 alpha2 walk: one_sheet and two_sheet both seen')

        call walk([1.0_rk, -1.0_rk, 1.0_rk, -1.4_rk, 0.0_rk], 5_ik, -0.6_rk, 0.01_rk, N4, 0.03_rk, &
                'R14 alpha3 walk', seen, inside, k4)
        call assert_true(.not. inside, 'R14 alpha3 walk: waist outside throughout')
    end subroutine test_r14_continuity_walks

    subroutine test_r15_r0_cross_checks()
        real(kind = rk), parameter :: TOL = 1.0e-5_rk
        type(geometry_t) :: geom

        call resolve_ok(ASYM_NECKED, geom, 'R15 ASYM_NECKED')
        call assert_abs_close(geom%u, 0.787101_rk, TOL, 'R15 ASYM u')
        call assert_abs_close(geom%z_join(1), -0.236690_rk, TOL, 'R15 ASYM z1')
        call assert_abs_close(geom%z_join(2), 0.557230_rk, TOL, 'R15 ASYM z2')
        call assert_abs_close(geom%l(3), 0.290020_rk, TOL, 'R15 ASYM l3')
        call assert_abs_close(geom%a_sq(3), 0.455200_rk, TOL, 'R15 ASYM a3^2')
        call assert_abs_close(geom%z_north, 1.500290_rk, TOL, 'R15 ASYM r_north')
        call assert_abs_close(-geom%z_south, 1.328790_rk, TOL, 'R15 ASYM r_south')
        call assert_abs_close(alpha1_f(geom), 0.134270_rk, TOL, 'R15 ASYM alpha1')

        call resolve_ok(PEAR, geom, 'R15 PEAR')
        call assert_abs_close(geom%u, 0.761898_rk, TOL, 'R15 PEAR u')
        call assert_abs_close(geom%z_join(1), -1.085240_rk, TOL, 'R15 PEAR z1')
        call assert_abs_close(geom%z_join(2), -0.704290_rk, TOL, 'R15 PEAR z2')
        call assert_abs_close(geom%z_north, 1.012910_rk, TOL, 'R15 PEAR r_north')
        call assert_abs_close(-geom%z_south, 1.159690_rk, TOL, 'R15 PEAR r_south')
        call assert_abs_close(alpha1_f(geom), -0.474390_rk, TOL, 'R15 PEAR alpha1')

        call resolve_ok(DEGENERATE_ONE, geom, 'R15 DEGENERATE_ONE')
        call assert_abs_close(geom%z_join(1), -0.894774_rk, TOL, 'R15 DEGENERATE_ONE z1')
        call assert_abs_close(geom%z_join(2), -0.894774_rk, TOL, 'R15 DEGENERATE_ONE z2')
    end subroutine test_r15_r0_cross_checks

    subroutine test_r16_symmetric_after_shift()
        real(kind = rk), parameter :: SYM(5, 3) = reshape([SYM_NECKED, INITIAL, SPHERE_B], [5, 3])
        type(geometry_t) :: geom
        integer(kind = ik) :: j
        character(len = 24) :: tag
        do j = 1_ik, 3_ik
            write(tag, '(A,I0)') 'R16 symmetric ', j
            call resolve_ok(SYM(:, j), geom, trim(tag))
            call assert_abs_close(geom%l(3), 0.0_rk, 1.0e-15_rk, trim(tag) // ': l3 = 0')
            call assert_abs_close(geom%z_south, -geom%z_north, 1.0e-15_rk, trim(tag) // ': z_south = -z_north')
        end do
    end subroutine test_r16_symmetric_after_shift

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
