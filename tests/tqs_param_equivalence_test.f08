!> Contract family 1, equivalence: the one-shot tier returns the same status
!! and the same bits as the cached tier for every output, on valid shapes
!! and on every rejection branch.
program tqs_param_equivalence_test
    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use tqs_parameterization_mod, only: cache_t, cache_radius_and_derivative_at_thetas_s, &
            compute_radius_and_derivative_standalone_s, TQS_ERROR_INVALID_COORDINATES, &
            TQS_ERROR_NO_GEOMETRY, TQS_ERROR_NOT_STAR_CONVEX, TQS_ERROR_SEPARATED, &
            SHAPE_ERROR_INVALID_GRID
    use tqs_test_fixtures_mod, only: SPHERE, SYM_NECKED, ASYM_NECKED, SPHEROID_MID, CYLINDER, &
            PEAR, SEPARATED, PINCHED, STAR_FAIL, BAD_K, NO_GEOM
    use tqs_test_tiers_mod, only: result_t, run_cached, run_one_shot, build_cache, &
            results_bits_eq, output_zero_filled, N_OUT, OUT_RADIUS_GRID
    use test_utils_mod, only: assert_true, assert_int_eq, arrays_bits_eq_f, all_zero_f, &
            test_summary
    implicit none

    integer(kind = ik), parameter :: NT = 32_ik, NZ = 257_ik, N_VALID = 9_ik
    real(kind = rk), parameter :: VALID_SET(5, N_VALID) = reshape([SPHERE, SYM_NECKED, &
            ASYM_NECKED, SPHEROID_MID, CYLINDER, PEAR, SEPARATED, PINCHED, STAR_FAIL], [5, N_VALID])
    character(len = 12), parameter :: VALID_NAMES(N_VALID) = [character(len = 12) :: &
            'SPHERE', 'SYM_NECKED', 'ASYM_NECKED', 'SPHEROID_MID', 'CYLINDER', 'PEAR', &
            'SEPARATED', 'PINCHED', 'STAR_FAIL']
    real(kind = rk) :: thetas(NT)
    type(cache_t) :: cache
    integer(kind = ik) :: i, status

    do i = 1_ik, NT
        thetas(i) = real(i, rk) * PI_C / real(NT + 1_ik, rk)
    end do
    call build_cache(thetas, cache, status)
    call assert_int_eq(status, 0_ik, 'cache built')

    call test_q0_negative_control()
    call test_q1_one_shot_equals_cached()
    call test_q2_rejections()

    call test_summary()

contains

    subroutine test_q0_negative_control()
        type(result_t) :: a, b
        call run_cached(cache, SYM_NECKED, thetas, NZ, a)
        call run_one_shot(ASYM_NECKED, thetas, NZ, b)
        call assert_true(.not. results_bits_eq(a, b), 'Q0 the comparator tells two shapes apart')
    end subroutine test_q0_negative_control

    subroutine test_q1_one_shot_equals_cached()
        type(result_t) :: a, b
        integer(kind = ik) :: f
        do f = 1_ik, N_VALID
            call run_cached(cache, VALID_SET(:, f), thetas, NZ, a)
            call run_one_shot(VALID_SET(:, f), thetas, NZ, b)
            call assert_true(results_bits_eq(a, b), 'Q1 ' // trim(VALID_NAMES(f)) // ': one-shot == cached, bitwise')
            call assert_true(arrays_bits_eq_f(a%radii_a, a%radii_d) .and. arrays_bits_eq_f(a%dr_a, a%dr_d), &
                    'Q1 ' // trim(VALID_NAMES(f)) // ': at_thetas == radius_and_derivative')
        end do
    end subroutine test_q1_one_shot_equals_cached

    !> Same code in both tiers on the rejected output, zero fill on both sides.
    subroutine rejection(params, slot, want, label)
        real(kind = rk), intent(in) :: params(5)
        integer(kind = ik), intent(in) :: slot, want
        character(len = *), intent(in) :: label
        type(result_t) :: a, b
        integer(kind = ik) :: j
        call run_cached(cache, params, thetas, NZ, a)
        call run_one_shot(params, thetas, NZ, b)
        call assert_int_eq(a%status(slot), want, label // ': cached code')
        call assert_int_eq(b%status(slot), want, label // ': one-shot code')
        call assert_true(output_zero_filled(a, slot) .and. output_zero_filled(b, slot), &
                label // ': zero-filled both sides')
        do j = 1_ik, N_OUT
            call assert_int_eq(a%status(j), b%status(j), label // ': every output, same status')
        end do
        call assert_true(results_bits_eq(a, b), label // ': every output, same bits')
    end subroutine rejection

    subroutine test_q2_rejections()
        real(kind = rk) :: bad_thetas(2), ra(2), da(2), rb(2), db(2)
        integer(kind = ik) :: sa, sb
        call rejection(BAD_K, OUT_RADIUS_GRID, TQS_ERROR_INVALID_COORDINATES, 'Q2 BAD_K')
        call rejection(NO_GEOM, OUT_RADIUS_GRID, TQS_ERROR_NO_GEOMETRY, 'Q2 NO_GEOM')
        call rejection(STAR_FAIL, OUT_RADIUS_GRID, TQS_ERROR_NOT_STAR_CONVEX, 'Q2 STAR_FAIL')
        call rejection(SEPARATED, OUT_RADIUS_GRID, TQS_ERROR_SEPARATED, 'Q2 SEPARATED')
        bad_thetas(1) = 0.5_rk
        bad_thetas(2) = 1.5_rk * PI_C
        call cache_radius_and_derivative_at_thetas_s(cache, SYM_NECKED, bad_thetas, ra, da, sa)
        call compute_radius_and_derivative_standalone_s(SYM_NECKED, bad_thetas, rb, db, sb)
        call assert_int_eq(sa, SHAPE_ERROR_INVALID_GRID, 'Q2 theta 3 pi / 2: cached code 3')
        call assert_int_eq(sb, SHAPE_ERROR_INVALID_GRID, 'Q2 theta 3 pi / 2: one-shot code 3')
        call assert_true(all_zero_f(ra) .and. all_zero_f(da) .and. all_zero_f(rb) .and. all_zero_f(db), &
                'Q2 theta 3 pi / 2: zero-filled both sides')
    end subroutine test_q2_rejections

end program tqs_param_equivalence_test
