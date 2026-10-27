!> Fixed-arity contract boundaries (L = TQS_N_PARAMS = 5): init arguments,
!! parameter counts in both tiers, per-position sensitivity, the profile
!! floor, output extents and the at-thetas theta count.
program tqs_param_boundary_test
    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use tqs_parameterization_mod, only: cache_t, tqs_geometry_t, cache_init_s, &
            cache_radius_grid_s, cache_radius_and_derivative_s, &
            cache_radius_and_derivative_at_thetas_s, cache_shape_s, cache_rho_z_grid_s, &
            cache_geometry_s, compute_radius_grid_standalone_s, &
            compute_radius_and_derivative_standalone_s, compute_shape_standalone_s, &
            compute_rho_z_grid_standalone_s, compute_geometry_standalone_s, &
            SHAPE_VALID, SHAPE_ERROR_TOO_MANY_PARAMS, SHAPE_ERROR_INVALID_GRID, &
            SHAPE_ERROR_WRONG_PARAM_COUNT, SHAPE_ERROR_INVALID_INIT
    use tqs_test_fixtures_mod, only: SYM_NECKED
    use test_utils_mod, only: assert_true, assert_int_eq, arrays_bits_eq_f, is_zero_f, &
            test_summary
    implicit none

    integer(kind = ik), parameter :: NT = 8_ik, NZ = 33_ik
    real(kind = rk), parameter :: NO_REALS(0) = [real(kind = rk) ::]
    real(kind = rk) :: thetas(NT)
    integer(kind = ik) :: i

    do i = 1_ik, NT
        thetas(i) = real(i, rk) * PI_C / real(NT + 1_ik, rk)
    end do

    call test_b1_init()
    call test_b2_cached_counts()
    call test_b3_one_shot_counts()
    call test_b4_every_position_matters()
    call test_b5_profile_floor()
    call test_b6_extents()
    call test_b7_at_thetas_count()

    call test_summary()

contains

    subroutine init_expect(max_params, init_thetas, want, label)
        integer(kind = ik), intent(in) :: max_params, want
        real(kind = rk), intent(in) :: init_thetas(:)
        character(len = *), intent(in) :: label
        type(cache_t) :: cache
        integer(kind = ik) :: status
        call cache_init_s(cache, max_params, init_thetas, status)
        call assert_int_eq(status, want, label)
    end subroutine init_expect

    subroutine test_b1_init()
        call init_expect(5_ik, thetas, SHAPE_VALID, 'B1 max_params 5 -> 0')
        call init_expect(6_ik, thetas, SHAPE_ERROR_TOO_MANY_PARAMS, 'B1 max_params 6 -> 1')
        call init_expect(4_ik, thetas, SHAPE_ERROR_INVALID_INIT, 'B1 max_params 4 -> 5')
        call init_expect(0_ik, thetas, SHAPE_ERROR_INVALID_INIT, 'B1 max_params 0 -> 5')
        call init_expect(5_ik, NO_REALS, SHAPE_ERROR_INVALID_GRID, 'B1 empty thetas -> 3')
        call init_expect(5_ik, [1.5_rk * PI_C], SHAPE_ERROR_INVALID_GRID, 'B1 theta 3 pi / 2 -> 3')
        call init_expect(5_ik, [PI_C], SHAPE_VALID, 'B1 theta pi -> 0')
        call init_expect(5_ik, [0.0_rk], SHAPE_VALID, 'B1 theta 0 -> 0')
    end subroutine test_b1_init

    !> Status of four cached outputs (radius_grid, shape, rho_z_grid, geometry).
    subroutine cached_counts(cache, params, want, label)
        type(cache_t), intent(in) :: cache
        real(kind = rk), intent(in) :: params(:)
        integer(kind = ik), intent(in) :: want
        character(len = *), intent(in) :: label
        real(kind = rk) :: radii(NT), rn, rs, z(NZ), rho(NZ), drho(NZ)
        type(tqs_geometry_t) :: geom
        integer(kind = ik) :: status
        call cache_radius_grid_s(cache, params, radii, status)
        call assert_int_eq(status, want, label // ' radius_grid')
        call cache_shape_s(cache, params, rn, rs, status)
        call assert_int_eq(status, want, label // ' shape')
        call cache_rho_z_grid_s(cache, params, z, rho, drho, status)
        call assert_int_eq(status, want, label // ' rho_z_grid')
        call cache_geometry_s(cache, params, geom, status)
        call assert_int_eq(status, want, label // ' geometry')
    end subroutine cached_counts

    subroutine one_shot_counts(params, want, label)
        real(kind = rk), intent(in) :: params(:)
        integer(kind = ik), intent(in) :: want
        character(len = *), intent(in) :: label
        real(kind = rk) :: radii(NT), rn, rs, z(NZ), rho(NZ), drho(NZ)
        type(tqs_geometry_t) :: geom
        integer(kind = ik) :: status
        call compute_radius_grid_standalone_s(params, thetas, radii, status)
        call assert_int_eq(status, want, label // ' radius_grid')
        call compute_shape_standalone_s(params, rn, rs, status)
        call assert_int_eq(status, want, label // ' shape')
        call compute_rho_z_grid_standalone_s(params, z, rho, drho, status)
        call assert_int_eq(status, want, label // ' rho_z_grid')
        call compute_geometry_standalone_s(params, geom, status)
        call assert_int_eq(status, want, label // ' geometry')
    end subroutine one_shot_counts

    subroutine test_b2_cached_counts()
        type(cache_t) :: cache
        integer(kind = ik) :: status
        real(kind = rk) :: p6(6)
        p6(1:5) = SYM_NECKED
        p6(6) = 0.0_rk
        call cache_init_s(cache, 5_ik, thetas, status)
        call cached_counts(cache, SYM_NECKED, SHAPE_VALID, 'B2 cached 5')
        call cached_counts(cache, p6, SHAPE_ERROR_WRONG_PARAM_COUNT, 'B2 cached 6')
        call cached_counts(cache, SYM_NECKED(1:4), SHAPE_ERROR_WRONG_PARAM_COUNT, 'B2 cached 4')
        call cached_counts(cache, NO_REALS, SHAPE_ERROR_WRONG_PARAM_COUNT, 'B2 cached 0')
    end subroutine test_b2_cached_counts

    subroutine test_b3_one_shot_counts()
        real(kind = rk) :: p6(6)
        p6(1:5) = SYM_NECKED
        p6(6) = 0.0_rk
        call one_shot_counts(SYM_NECKED, SHAPE_VALID, 'B3 one-shot 5')
        call one_shot_counts(p6, SHAPE_ERROR_TOO_MANY_PARAMS, 'B3 one-shot 6')
        call one_shot_counts(SYM_NECKED(1:4), SHAPE_ERROR_WRONG_PARAM_COUNT, 'B3 one-shot 4')
        call one_shot_counts(NO_REALS, SHAPE_ERROR_WRONG_PARAM_COUNT, 'B3 one-shot 0')
    end subroutine test_b3_one_shot_counts

    subroutine test_b4_every_position_matters()
        type(cache_t) :: cache
        real(kind = rk) :: base(NT), moved(NT), p(5)
        integer(kind = ik) :: status, slot
        character(len = 24) :: tag
        call cache_init_s(cache, 5_ik, thetas, status)
        call cache_radius_grid_s(cache, SYM_NECKED, base, status)
        call assert_int_eq(status, SHAPE_VALID, 'B4 base status 0')
        do slot = 1_ik, 5_ik
            write(tag, '(A,I0)') 'B4 position ', slot
            p = SYM_NECKED
            p(slot) = p(slot) + 0.05_rk
            call cache_radius_grid_s(cache, p, moved, status)
            call assert_int_eq(status, SHAPE_VALID, trim(tag) // ': status 0')
            call assert_true(.not. arrays_bits_eq_f(base, moved), trim(tag) // ': changes radius_grid')
        end do
    end subroutine test_b4_every_position_matters

    subroutine test_b5_profile_floor()
        type(cache_t) :: cache
        real(kind = rk) :: z3(3), r3(3), d3(3), z2(2), r2(2), d2(2)
        integer(kind = ik) :: status
        call cache_init_s(cache, 5_ik, thetas, status)
        call cache_rho_z_grid_s(cache, SYM_NECKED, z3, r3, d3, status)
        call assert_int_eq(status, SHAPE_VALID, 'B5 cached size 3 -> 0')
        call assert_true(is_zero_f(r3(1)) .and. is_zero_f(r3(3)) .and. r3(2) > 0.0_rk, &
                'B5 cached size 3: poles 0, middle positive')
        call cache_rho_z_grid_s(cache, SYM_NECKED, z2, r2, d2, status)
        call assert_int_eq(status, SHAPE_ERROR_INVALID_GRID, 'B5 cached size 2 -> 3')
        call compute_rho_z_grid_standalone_s(SYM_NECKED, z3, r3, d3, status)
        call assert_int_eq(status, SHAPE_VALID, 'B5 one-shot size 3 -> 0')
        call compute_rho_z_grid_standalone_s(SYM_NECKED, z2, r2, d2, status)
        call assert_int_eq(status, SHAPE_ERROR_INVALID_GRID, 'B5 one-shot size 2 -> 3')
    end subroutine test_b5_profile_floor

    subroutine test_b6_extents()
        type(cache_t) :: cache
        real(kind = rk) :: radii(NT + 1_ik), dr(NT + 1_ik), z(NZ), rho(NZ), drho(NZ - 1_ik)
        integer(kind = ik) :: status
        call cache_init_s(cache, 5_ik, thetas, status)
        call cache_radius_grid_s(cache, SYM_NECKED, radii, status)
        call assert_int_eq(status, SHAPE_ERROR_INVALID_GRID, 'B6 radius_grid n_thetas + 1 -> 3')
        call cache_radius_and_derivative_s(cache, SYM_NECKED, radii, dr, status)
        call assert_int_eq(status, SHAPE_ERROR_INVALID_GRID, 'B6 radius_and_derivative n_thetas + 1 -> 3')
        call cache_radius_and_derivative_at_thetas_s(cache, SYM_NECKED, thetas, radii, dr, status)
        call assert_int_eq(status, SHAPE_ERROR_INVALID_GRID, 'B6 at_thetas size(radii) /= size(thetas) -> 3')
        call cache_rho_z_grid_s(cache, SYM_NECKED, z, rho, drho, status)
        call assert_int_eq(status, SHAPE_ERROR_INVALID_GRID, 'B6 size(drho_dz) /= size(z) -> 3')
        call compute_radius_grid_standalone_s(SYM_NECKED, thetas, radii, status)
        call assert_int_eq(status, SHAPE_ERROR_INVALID_GRID, 'B6 one-shot radius_grid extent -> 3')
        call compute_radius_and_derivative_standalone_s(SYM_NECKED, thetas, radii, dr, status)
        call assert_int_eq(status, SHAPE_ERROR_INVALID_GRID, 'B6 one-shot radius_and_derivative extent -> 3')
    end subroutine test_b6_extents

    subroutine test_b7_at_thetas_count()
        type(cache_t) :: cache
        real(kind = rk) :: r1(1), d1(1), r0(0), d0(0)
        integer(kind = ik) :: status
        call cache_init_s(cache, 5_ik, thetas, status)
        call cache_radius_and_derivative_at_thetas_s(cache, SYM_NECKED, [0.7_rk], r1, d1, status)
        call assert_int_eq(status, SHAPE_VALID, 'B7 at_thetas one theta -> 0')
        call cache_radius_and_derivative_at_thetas_s(cache, SYM_NECKED, NO_REALS, r0, d0, status)
        call assert_int_eq(status, SHAPE_ERROR_INVALID_GRID, 'B7 at_thetas no theta -> 3')
    end subroutine test_b7_at_thetas_count

end program tqs_param_boundary_test
