!> Contract family 2, statelessness: a call's result depends only on its
!! arguments — not on earlier calls on the same cache, not on which cache,
!! and not on the caller's theta array after init.
program tqs_param_statelessness_test
    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use tqs_parameterization_mod, only: cache_t
    use tqs_test_fixtures_mod, only: SYM_NECKED, ASYM_NECKED, SEPARATED, BAD_K, NO_GEOM
    use tqs_test_tiers_mod, only: result_t, run_cached, build_cache, results_bits_eq
    use test_utils_mod, only: assert_true, assert_int_eq, test_summary
    implicit none

    integer(kind = ik), parameter :: NT = 32_ik, NZ = 257_ik
    real(kind = rk) :: thetas(NT)
    integer(kind = ik) :: i

    do i = 1_ik, NT
        thetas(i) = real(i, rk) * PI_C / real(NT + 1_ik, rk)
    end do

    call test_s0_negative_control()
    call test_s1_repeat()
    call test_s2_after_rejection()
    call test_s3_cache_identity()

    call test_summary()

contains

    subroutine fresh(cache)
        type(cache_t), intent(out) :: cache
        integer(kind = ik) :: status
        call build_cache(thetas, cache, status)
        call assert_int_eq(status, 0_ik, 'cache built')
    end subroutine fresh

    subroutine test_s0_negative_control()
        type(cache_t) :: cache
        type(result_t) :: a, b
        call fresh(cache)
        call run_cached(cache, SYM_NECKED, thetas, NZ, a)
        call run_cached(cache, ASYM_NECKED, thetas, NZ, b)
        call assert_true(.not. results_bits_eq(a, b), 'S0 the comparator tells two shapes apart')
    end subroutine test_s0_negative_control

    subroutine test_s1_repeat()
        type(cache_t) :: cache
        type(result_t) :: a, b
        call fresh(cache)
        call run_cached(cache, ASYM_NECKED, thetas, NZ, a)
        call run_cached(cache, ASYM_NECKED, thetas, NZ, b)
        call assert_true(results_bits_eq(a, b), 'S1 repeat on one cache: identical bits')
    end subroutine test_s1_repeat

    subroutine after(between, label)
        real(kind = rk), intent(in) :: between(:)
        character(len = *), intent(in) :: label
        type(cache_t) :: used, clean
        type(result_t) :: reference, sink, got
        call fresh(clean)
        call run_cached(clean, ASYM_NECKED, thetas, NZ, reference)
        call fresh(used)
        call run_cached(used, between, thetas, NZ, sink)
        call run_cached(used, ASYM_NECKED, thetas, NZ, got)
        call assert_true(results_bits_eq(reference, got), 'S2 after ' // label // ': same bits as a fresh cache')
    end subroutine after

    subroutine test_s2_after_rejection()
        real(kind = rk) :: p4(4)
        p4 = ASYM_NECKED(1:4)
        call after(BAD_K, 'BAD_K')
        call after(NO_GEOM, 'NO_GEOM')
        call after(SEPARATED, 'SEPARATED')
        call after(p4, 'a 4-element vector')
    end subroutine test_s2_after_rejection

    subroutine test_s3_cache_identity()
        type(cache_t) :: c1, c2
        type(result_t) :: a, b, c
        real(kind = rk) :: own(NT), at_thetas(NT)
        integer(kind = ik) :: status
        call fresh(c1)
        call fresh(c2)
        call run_cached(c1, ASYM_NECKED, thetas, NZ, a)
        call run_cached(c2, ASYM_NECKED, thetas, NZ, b)
        call assert_true(results_bits_eq(a, b), 'S3 two caches on the same thetas: identical bits')

        ! The cache owns a copy: overwrite the caller's array after init.
        own = thetas
        at_thetas = thetas
        call build_cache(own, c1, status)
        call assert_int_eq(status, 0_ik, 'S3 cache on the caller array built')
        own = 0.25_rk
        call run_cached(c1, ASYM_NECKED, at_thetas, NZ, c)
        call assert_true(results_bits_eq(a, c), 'S3 caller thetas overwritten after init: results do not move')
    end subroutine test_s3_cache_identity

end program tqs_param_statelessness_test
