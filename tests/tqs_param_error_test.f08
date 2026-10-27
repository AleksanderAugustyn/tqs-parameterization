!> The full status matrix of the ten outputs in both tiers: valid rows,
!! value codes 100-103 with zero fill, usage-code precedence, and the
!! identities between outputs (at-thetas == grid, unchecked == checked).
program tqs_param_error_test
    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use tqs_parameterization_mod, only: cache_t, cache_free_s, &
            cache_radius_and_derivative_at_thetas_s, cache_rho_z_grid_s, &
            compute_radius_grid_standalone_s, compute_rho_z_grid_standalone_s, &
            TQS_ERROR_INVALID_COORDINATES, TQS_ERROR_NO_GEOMETRY, TQS_ERROR_NOT_STAR_CONVEX, &
            TQS_ERROR_SEPARATED, SHAPE_ERROR_WRONG_PARAM_COUNT, SHAPE_ERROR_TOO_MANY_PARAMS, &
            SHAPE_ERROR_INVALID_GRID, SHAPE_ERROR_CACHE_NOT_INITIALIZED
    use tqs_test_fixtures_mod, only: SPHERE, SPHERE_B, SYM_NECKED, ASYM_NECKED, PEAR, PINCHED, &
            SEPARATED, STAR_FAIL, BAD_K, NO_GEOM
    use tqs_test_tiers_mod, only: result_t, run_cached, run_one_shot, build_cache, &
            output_zero_filled, geometries_bits_eq, N_OUT, OUT_RADIUS_GRID, OUT_SHAPE, &
            OUT_RHO_Z, OUT_SCISSION
    use test_utils_mod, only: assert_true, assert_int_eq, assert_abs_close, bits_eq_f, &
            arrays_bits_eq_f, is_zero_f, all_zero_f, test_summary
    implicit none

    integer(kind = ik), parameter :: NT = 8_ik, NZ = 101_ik, MID = 51_ik
    character(len = 9), parameter :: TIER_NAMES(2) = [character(len = 9) :: 'cached', 'one-shot']
    real(kind = rk) :: thetas(NT)
    integer(kind = ik) :: i

    do i = 1_ik, NT
        thetas(i) = real(i, rk) * PI_C / real(NT + 1_ik, rk)
    end do

    call test_e1_valid()
    call test_e2_value_codes()
    call test_e3_separated()
    call test_e4_pinched()
    call test_e5_star_fail()
    call test_e6_precedence()
    call test_e7_at_thetas_identity()
    call test_e8_identities()

    call test_summary()

contains

    !> Both tiers of one parameter vector.
    subroutine run_both(params, res)
        real(kind = rk), intent(in) :: params(:)
        type(result_t), intent(out) :: res(2)
        type(cache_t) :: cache
        integer(kind = ik) :: status
        call build_cache(thetas, cache, status)
        call assert_int_eq(status, 0_ik, 'cache built')
        call run_cached(cache, params, thetas, NZ, res(1))
        call run_one_shot(params, thetas, NZ, res(2))
    end subroutine run_both

    !> Statuses as expected, and every output with a nonzero status zero-filled.
    subroutine expect_statuses(res, want, label)
        type(result_t), intent(in) :: res
        integer(kind = ik), intent(in) :: want(N_OUT)
        character(len = *), intent(in) :: label
        integer(kind = ik) :: j
        character(len = 8) :: tag
        do j = 1_ik, N_OUT
            write(tag, '(A,I0)') ' out ', j
            call assert_int_eq(res%status(j), want(j), label // tag // ' status')
            if (want(j) /= 0_ik) call assert_true(output_zero_filled(res, j), &
                    label // tag // ' zero-filled')
        end do
    end subroutine expect_statuses

    subroutine test_e1_valid()
        type(result_t) :: res(2)
        real(kind = rk) :: poles(2)
        integer(kind = ik) :: t, status, want(N_OUT)
        call run_both(SYM_NECKED, res)
        call compute_radius_grid_standalone_s(SYM_NECKED, [0.0_rk, PI_C], poles, status)
        call assert_int_eq(status, 0_ik, 'E1 SYM_NECKED pole rays status 0')
        want = 0_ik
        do t = 1_ik, 2_ik
            call expect_statuses(res(t), want, 'E1 SYM_NECKED ' // trim(TIER_NAMES(t)))
            call assert_true(bits_eq_f(res(t)%r_north, res(t)%geom%r_north) .and. &
                    bits_eq_f(res(t)%r_south, res(t)%geom%r_south), &
                    'E1 ' // trim(TIER_NAMES(t)) // ': shape poles = geometry poles')
            call assert_abs_close(res(t)%r_north, poles(1), 1.0e-14_rk, 'E1 ' // trim(TIER_NAMES(t)) // ': r_north = R(0)')
            call assert_abs_close(res(t)%r_south, poles(2), 1.0e-14_rk, 'E1 ' // trim(TIER_NAMES(t)) // ': r_south = R(pi)')
        end do
    end subroutine test_e1_valid

    subroutine test_e2_value_codes()
        type(result_t) :: res(2)
        integer(kind = ik) :: t, want(N_OUT)
        call run_both(BAD_K, res)
        want = TQS_ERROR_INVALID_COORDINATES
        do t = 1_ik, 2_ik
            call expect_statuses(res(t), want, 'E2 BAD_K ' // trim(TIER_NAMES(t)))
        end do
        call run_both(NO_GEOM, res)
        want = TQS_ERROR_NO_GEOMETRY
        want(OUT_SCISSION) = 0_ik
        do t = 1_ik, 2_ik
            call expect_statuses(res(t), want, 'E2 NO_GEOM ' // trim(TIER_NAMES(t)))
            call assert_true(.not. res(t)%defined .and. res(t)%s1t > 0.0_rk, &
                    'E2 NO_GEOM ' // trim(TIER_NAMES(t)) // ': scission valid, sigma1_star undefined')
        end do
    end subroutine test_e2_value_codes

    subroutine test_e3_separated()
        type(result_t) :: res(2)
        integer(kind = ik) :: t, want(N_OUT)
        character(len = 24) :: tag
        call run_both(SEPARATED, res)
        want = 0_ik
        want(OUT_RADIUS_GRID:OUT_RHO_Z) = TQS_ERROR_SEPARATED
        do t = 1_ik, 2_ik
            tag = 'E3 SEPARATED ' // trim(TIER_NAMES(t))
            call expect_statuses(res(t), want, trim(tag))
            call assert_true(is_zero_f(res(t)%rhou(MID)) .and. res(t)%rhou(2) > 0.0_rk &
                    .and. res(t)%rhou(NZ - 1_ik) > 0.0_rk, trim(tag) // ': unchecked draws both fragments')
            call assert_true(res(t)%found .and. is_zero_f(res(t)%rho_neck), trim(tag) // ': neck at the split')
            call assert_true(res(t)%geom%is_separated .and. res(t)%geom%r_north > 0.0_rk &
                    .and. res(t)%geom%r_south > 0.0_rk, trim(tag) // ': geometry flags and poles')
            call assert_true(is_zero_f(res(t)%d), trim(tag) // ': d = 0')
            call assert_true(res(t)%defined, trim(tag) // ': scission defined')
        end do
    end subroutine test_e3_separated

    subroutine test_e4_pinched()
        type(result_t) :: res(2)
        integer(kind = ik) :: t, want(N_OUT)
        character(len = 24) :: tag
        call run_both(PINCHED, res)
        want = 0_ik
        want(OUT_RADIUS_GRID:OUT_SHAPE) = TQS_ERROR_NOT_STAR_CONVEX
        do t = 1_ik, 2_ik
            tag = 'E4 PINCHED ' // trim(TIER_NAMES(t))
            call expect_statuses(res(t), want, trim(tag))
            call assert_true(arrays_bits_eq_f(res(t)%z, res(t)%zu) .and. arrays_bits_eq_f(res(t)%rho, res(t)%rhou) &
                    .and. arrays_bits_eq_f(res(t)%drho, res(t)%drhou), trim(tag) // ': checked = unchecked')
            call assert_true(res(t)%found .and. is_zero_f(res(t)%rho_neck), trim(tag) // ': neck rho = 0')
            call assert_true(res(t)%geom%is_pinched, trim(tag) // ': geometry is_pinched')
            call assert_true(is_zero_f(res(t)%d), trim(tag) // ': d = 0')
        end do
    end subroutine test_e4_pinched

    subroutine test_e5_star_fail()
        type(result_t) :: res(2)
        integer(kind = ik) :: t, want(N_OUT)
        character(len = 24) :: tag
        call run_both(STAR_FAIL, res)
        want = 0_ik
        want(OUT_RADIUS_GRID:OUT_SHAPE) = TQS_ERROR_NOT_STAR_CONVEX
        do t = 1_ik, 2_ik
            tag = 'E5 STAR_FAIL ' // trim(TIER_NAMES(t))
            call expect_statuses(res(t), want, trim(tag))
            call assert_true(res(t)%found .and. res(t)%rho_neck > 0.0_rk, trim(tag) // ': neck found')
            call assert_true(.not. res(t)%geom%is_separated .and. .not. res(t)%geom%is_pinched, &
                    trim(tag) // ': connected, not pinched')
        end do
    end subroutine test_e5_star_fail

    subroutine test_e6_precedence()
        real(kind = rk) :: p4(4), p6(6), radii(NT), dr(NT), z(NZ), rho(NZ - 1_ik), drho(NZ)
        type(result_t) :: res
        type(cache_t) :: cache
        integer(kind = ik) :: status, want(N_OUT)

        p4 = SEPARATED(1:4)
        p6(1:5) = SEPARATED
        p6(6) = 0.0_rk
        call build_cache(thetas, cache, status)
        want = SHAPE_ERROR_WRONG_PARAM_COUNT
        call run_cached(cache, p4, thetas, NZ, res)
        call expect_statuses(res, want, 'E6 SEPARATED 4 params cached')
        call run_one_shot(p4, thetas, NZ, res)
        call expect_statuses(res, want, 'E6 SEPARATED 4 params one-shot')
        call run_cached(cache, p6, thetas, NZ, res)
        call expect_statuses(res, want, 'E6 SEPARATED 6 params cached')
        want = SHAPE_ERROR_TOO_MANY_PARAMS
        call run_one_shot(p6, thetas, NZ, res)
        call expect_statuses(res, want, 'E6 SEPARATED 6 params one-shot')

        call cache_radius_and_derivative_at_thetas_s(cache, PINCHED, [0.5_rk, 1.5_rk * PI_C], &
                radii(1:2), dr(1:2), status)
        call assert_int_eq(status, SHAPE_ERROR_INVALID_GRID, 'E6 PINCHED at_thetas with 3 pi / 2 -> 3')
        call assert_true(all_zero_f(radii(1:2)) .and. all_zero_f(dr(1:2)), 'E6 bad theta zero-filled')

        call cache_rho_z_grid_s(cache, SEPARATED, z, rho, drho, status)
        call assert_int_eq(status, SHAPE_ERROR_INVALID_GRID, 'E6 SEPARATED cached rho_z unequal sizes -> 3')
        call compute_rho_z_grid_standalone_s(SEPARATED, z, rho, drho, status)
        call assert_int_eq(status, SHAPE_ERROR_INVALID_GRID, 'E6 SEPARATED one-shot rho_z unequal sizes -> 3')

        call cache_free_s(cache)
        want = SHAPE_ERROR_CACHE_NOT_INITIALIZED
        call run_cached(cache, BAD_K, thetas, NZ, res)
        call expect_statuses(res, want, 'E6 freed cache with BAD_K -> 2')
    end subroutine test_e6_precedence

    subroutine test_e7_at_thetas_identity()
        type(result_t) :: res
        type(cache_t) :: cache
        integer(kind = ik) :: status
        call build_cache(thetas, cache, status)
        call run_cached(cache, SYM_NECKED, thetas, NZ, res)
        call assert_true(arrays_bits_eq_f(res%radii_a, res%radii_d) .and. arrays_bits_eq_f(res%dr_a, res%dr_d), &
                'E7 SYM_NECKED at_thetas = radius_and_derivative bitwise')
        call assert_true(arrays_bits_eq_f(res%radii_g, res%radii_d), 'E7 SYM_NECKED radius_grid = radius_and_derivative radii')
        call run_cached(cache, PEAR, thetas, NZ, res)
        call assert_true(arrays_bits_eq_f(res%radii_a, res%radii_d) .and. arrays_bits_eq_f(res%dr_a, res%dr_d), &
                'E7 PEAR at_thetas = radius_and_derivative bitwise')
    end subroutine test_e7_at_thetas_identity

    subroutine test_e8_identities()
        type(result_t) :: res(2)
        integer(kind = ik) :: t
        call run_both(ASYM_NECKED, res)
        do t = 1_ik, 2_ik
            call assert_true(arrays_bits_eq_f(res(t)%z, res(t)%zu) .and. arrays_bits_eq_f(res(t)%rho, res(t)%rhou) &
                    .and. arrays_bits_eq_f(res(t)%drho, res(t)%drhou), &
                    'E8 ASYM_NECKED ' // trim(TIER_NAMES(t)) // ': unchecked = checked bitwise')
        end do
        call run_both(SPHERE, res)
        call assert_true(maxval(abs(res(1)%radii_g - 1.0_rk)) <= 1.0e-14_rk, 'E8 SPHERE radius_grid = 1')
        call run_both(SPHERE_B, res)
        call assert_true(maxval(abs(res(1)%radii_g - 1.0_rk)) <= 1.0e-14_rk, 'E8 SPHERE_B radius_grid = 1')
        call assert_true(geometries_bits_eq(res(1)%geom, res(2)%geom), 'E8 SPHERE_B geometry tiers agree')
    end subroutine test_e8_identities

end program tqs_param_error_test
