!> Breakdown: where the closed form needs a rule or loses accuracy. The
!! degenerate-curvature tolerances, the sigma2 = 0 and sigma1 = 0 crossings,
!! |alpha2| -> 2, a3^2 through 0 with the waist inside and outside, large
!! coordinates, the star-margin edge and the both-pass root count. Each block
!! prints what it records; the assertions are the plan's.
program tqs_param_breakdown_test
    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use tqs_workers_mod, only: geometry_t, body_integrals_s, star_margin_f
    use tqs_resolve_mod, only: resolve_unit_s, resolve_s
    use tqs_parameterization_mod, only: compute_radius_grid_standalone_s, &
            compute_neck_standalone_s, TQS_DEGENERATE_TOL, STAR_CONVEXITY_MARGIN, &
            TQS_MID_SPHEROID, TQS_MID_CYLINDER, TQS_MID_ONE_SHEET, &
            TQS_ERROR_INVALID_COORDINATES, TQS_ERROR_NO_GEOMETRY, &
            TQS_ERROR_NOT_STAR_CONVEX, TQS_ERROR_SEPARATED
    use tqs_test_fixtures_mod, only: SPHERE_B, pv
    use tqs_test_tiers_mod, only: result_t, run_one_shot, N_OUT
    use test_utils_mod, only: assert_true, assert_int_eq, bits_eq_f, test_summary
    implicit none

    integer(kind = ik), parameter :: NT = 16_ik
    real(kind = rk) :: thetas(NT)
    integer(kind = ik) :: i

    do i = 1_ik, NT
        thetas(i) = real(i, rk) * PI_C / real(NT + 1_ik, rk)
    end do

    call test_k1_coincident_family()
    call test_k2_non_coincident()
    call test_k3_both_ends()
    call test_k4_sigma2_to_zero()
    call test_k5_sigma1_to_zero()
    call test_k6_alpha2_to_two()
    call test_k7_a3_through_zero()
    call test_k8_sphere_has_no_neck()
    call test_k9_large_coordinates()
    call test_k10_star_margin_edge()
    call test_k11_both_pass_count()

    call test_summary()

contains

    subroutine test_k1_coincident_family()
        real(kind = rk) :: p(5)
        integer(kind = ik), parameter :: N = 15_ik
        real(kind = rk), parameter :: DS(N) = [1.0e-1_rk, 1.0e-2_rk, 1.0e-3_rk, 1.0e-4_rk, &
                1.0e-5_rk, 1.0e-6_rk, 1.0e-8_rk, 1.0e-10_rk, 1.0e-12_rk, 0.0_rk, -1.0e-5_rk, &
                -1.0e-4_rk, -1.0e-3_rk, -1.0e-2_rk, -1.0e-1_rk]
        type(geometry_t) :: g, g0
        integer(kind = ik) :: j, status, n_valid
        character(len = 40) :: tag
        p = pv(0.5_rk, 2.0_rk, 1.5_rk, -0.5_rk, 1.0_rk)
        call resolve_s(p, g0, status)
        write(*, '(A)') 'K1 coincident family [0.5, 2 + d, 1.5, -0.5, 1]: d, status, n_valid, z1, z2'
        do j = 1_ik, N
            write(tag, '(A,ES9.1)') 'K1 d = ', DS(j)
            p = pv(0.5_rk, 2.0_rk + DS(j), 1.5_rk, -0.5_rk, 1.0_rk)
            call resolve_s(p, g, status, n_valid)
            write(*, '(A,ES10.2,2I5,2ES24.16)') '   ', DS(j), status, n_valid, g%z_join(1), g%z_join(2)
            call assert_int_eq(status, 0_ik, trim(tag) // ': valid')
            call assert_true(abs(g%z_join(1) - g0%z_join(1)) <= 1.0e-4_rk, trim(tag) // ': z1 continuous')
            if (abs(DS(j)) < TQS_DEGENERATE_TOL) then
                call assert_true(bits_eq_f(g%z_join(1), g%z_join(2)), trim(tag) // ': degenerate path, z1 = z2 bitwise')
                call assert_int_eq(n_valid, 0_ik, trim(tag) // ': degenerate path, n_valid_roots = 0')
            end if
            if (DS(j) > 5.0e-3_rk) call assert_int_eq(n_valid, 1_ik, trim(tag) // ': generic path, one valid root')
        end do
    end subroutine test_k1_coincident_family

    subroutine test_k2_non_coincident()
        real(kind = rk) :: p(5)
        real(kind = rk), parameter :: DS(5) = [1.0e-1_rk, 1.0e-2_rk, 1.0e-4_rk, 1.0e-8_rk, 0.0_rk]
        type(geometry_t) :: g
        integer(kind = ik) :: j, status
        character(len = 40) :: tag
        do j = 1_ik, 5_ik
            write(tag, '(A,ES9.1)') 'K2 d = ', DS(j)
            p = pv(0.5_rk, 2.0_rk + DS(j), 1.5_rk, -0.3_rk, 1.0_rk)
            call resolve_s(p, g, status)
            call assert_int_eq(status, TQS_ERROR_NO_GEOMETRY, trim(tag) // ': non-coincident -> 101')
        end do
    end subroutine test_k2_non_coincident

    subroutine test_k3_both_ends()
        real(kind = rk) :: p(5)
        real(kind = rk), parameter :: DS(3) = [0.0_rk, 1.0e-5_rk, -1.0e-5_rk]
        type(geometry_t) :: g
        integer(kind = ik) :: j, status
        character(len = 40) :: tag
        do j = 1_ik, 3_ik
            write(tag, '(A,ES9.1)') 'K3 d = ', DS(j)
            p = pv(0.0_rk, 1.5_rk + DS(j), 1.5_rk, 0.0_rk, 0.0_rk)
            call resolve_s(p, g, status)
            call assert_int_eq(status, 0_ik, trim(tag) // ': spheroid line valid')
        end do
        p = pv(1.0e-3_rk, 1.5_rk, 1.5_rk, 0.0_rk, 0.0_rk)
        call resolve_s(p, g, status)
        call assert_int_eq(status, TQS_ERROR_NO_GEOMETRY, 'K3 sigma1 = 1e-3 on the line -> 101')
        p = pv(0.0_rk, 1.5_rk, 1.5_rk, 1.0e-3_rk, 0.0_rk)
        call resolve_s(p, g, status)
        call assert_int_eq(status, TQS_ERROR_NO_GEOMETRY, 'K3 alpha2 = 1e-3 on the line -> 101')
    end subroutine test_k3_both_ends

    subroutine test_k4_sigma2_to_zero()
        real(kind = rk) :: p(5)
        real(kind = rk), parameter :: S2(7) = [1.0e-2_rk, -1.0e-2_rk, 1.0e-4_rk, -1.0e-4_rk, &
                1.0e-8_rk, -1.0e-8_rk, 0.0_rk]
        type(geometry_t) :: g, g0
        integer(kind = ik) :: j, status
        real(kind = rk) :: tol
        character(len = 40) :: tag
        p = pv(0.6_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.3_rk)
        call resolve_s(p, g0, status)
        call assert_int_eq(g0%mid_kind, TQS_MID_CYLINDER, 'K4 sigma2 = 0: cylinder')
        write(*, '(A)') 'K4 [0.6, s2, 1, 0, 0.3]: s2, kind, l3, a3^2'
        do j = 1_ik, 7_ik
            write(tag, '(A,ES9.1)') 'K4 s2 = ', S2(j)
            p = pv(0.6_rk, S2(j), 1.0_rk, 0.0_rk, 0.3_rk)
            call resolve_s(p, g, status)
            write(*, '(A,ES10.2,I3,2ES24.16)') '   ', S2(j), g%mid_kind, g%l(3), g%a_sq(3)
            call assert_int_eq(status, 0_ik, trim(tag) // ': valid')
            if (S2(j) > 0.0_rk) call assert_int_eq(g%mid_kind, TQS_MID_SPHEROID, trim(tag) // ': spheroid')
            if (S2(j) < 0.0_rk) call assert_int_eq(g%mid_kind, TQS_MID_ONE_SHEET, trim(tag) // ': one_sheet')
            ! l3 and a3^2 move at first order in sigma2 (slopes ~0.023 and ~0.09):
            ! 1e-3 at |sigma2| = 1e-2 as the plan, scaled linearly below it.
            tol = 0.1_rk * abs(S2(j))
            call assert_true(abs(g%l(3) - g0%l(3)) <= tol .and. abs(g%a_sq(3) - g0%a_sq(3)) <= tol, &
                    trim(tag) // ': l3 and a3^2 continuous through 0')
        end do

        write(*, '(A)') 'K4 [0.6, s2, 1, 0.2, 0.3]: s2, status, l3, a3^2 (alpha2 /= 0: l3, a3^2 diverge at 0)'
        do j = 1_ik, 2_ik
            p = pv(0.6_rk, S2(j), 1.0_rk, 0.2_rk, 0.3_rk)
            call resolve_s(p, g, status)
            write(*, '(A,ES10.2,I5,2ES24.16)') '   ', S2(j), status, g%l(3), g%a_sq(3)
            call assert_int_eq(status, 0_ik, 'K4 alpha2 = 0.2 valid at sigma2 = +-1e-2')
        end do
        p = pv(0.6_rk, 0.0_rk, 1.0_rk, 0.2_rk, 0.3_rk)
        call resolve_s(p, g, status)
        call assert_int_eq(status, TQS_ERROR_NO_GEOMETRY, 'K4 alpha2 = 0.2, sigma2 = 0 -> 101')
    end subroutine test_k4_sigma2_to_zero

    subroutine test_k5_sigma1_to_zero()
        real(kind = rk) :: p(5)
        real(kind = rk), parameter :: S1(6) = [1.0e-1_rk, 1.0e-2_rk, 1.0e-4_rk, 1.0e-8_rk, 1.0e-12_rk, 0.0_rk]
        real(kind = rk), parameter :: S1_BAD(3) = [0.1_rk, 0.01_rk, 0.0_rk]
        type(geometry_t) :: g, gu
        integer(kind = ik) :: j, status
        real(kind = rk) :: prev
        character(len = 40) :: tag
        prev = huge(1.0_rk)
        do j = 1_ik, 6_ik
            write(tag, '(A,ES9.1)') 'K5 s1 = ', S1(j)
            p = pv(S1(j), -0.6_rk, 1.0_rk, 0.0_rk, 0.3_rk)
            call resolve_s(p, g, status)
            call assert_int_eq(status, 0_ik, trim(tag) // ': valid')
            call assert_true(g%u <= 1.2_rk .and. g%z_join(2) - g%z_join(1) <= 2.0_rk * S1(j), &
                    trim(tag) // ': middle no longer than 2 sigma1')
            p = pv(S1(j), -0.6_rk, 1.0_rk, 0.0_rk, 0.3_rk)
            call resolve_unit_s(p, gu, status)
            call assert_true(max(abs(gu%z_join(1)), abs(gu%z_join(2))) <= prev, &
                    trim(tag) // ': unit-stage joins shrink monotonically')
            prev = max(abs(gu%z_join(1)), abs(gu%z_join(2)))
        end do
        call assert_true(prev <= 0.0_rk, 'K5 sigma1 = 0: unit-stage joins at 0')
        do j = 1_ik, 3_ik
            write(tag, '(A,ES9.1)') 'K5 alpha2 = 0.3, s1 = ', S1_BAD(j)
            p = pv(S1_BAD(j), -0.6_rk, 1.0_rk, 0.3_rk, 0.0_rk)
            call resolve_s(p, g, status)
            call assert_int_eq(status, TQS_ERROR_NO_GEOMETRY, trim(tag) // ': unbridgeable -> 101')
        end do
    end subroutine test_k5_sigma1_to_zero

    subroutine test_k6_alpha2_to_two()
        real(kind = rk) :: p(5)
        real(kind = rk), parameter :: A2(4) = [1.9_rk, 1.99_rk, 1.999_rk, 1.9999_rk]
        type(geometry_t) :: g
        integer(kind = ik) :: j, status
        write(*, '(A)') 'K6 [0.5, -0.6, 1, a2, 0]: a2, status'
        do j = 1_ik, 4_ik
            p = pv(0.5_rk, -0.6_rk, 1.0_rk, A2(j), 0.0_rk)
            call resolve_s(p, g, status)
            write(*, '(A,F8.4,I5)') '   ', A2(j), status
            call assert_true(status == 0_ik .or. status == TQS_ERROR_NO_GEOMETRY, &
                    'K6 |alpha2| -> 2: valid or 101, never another code')
        end do
        p = pv(0.5_rk, -0.6_rk, 1.0_rk, 2.0_rk, 0.0_rk)
        call resolve_s(p, g, status)
        call assert_int_eq(status, TQS_ERROR_INVALID_COORDINATES, 'K6 alpha2 = 2 -> 100')
    end subroutine test_k6_alpha2_to_two

    subroutine test_k7_a3_through_zero()
        real(kind = rk) :: p(5)
        real(kind = rk), parameter :: ES(3) = [-1.0e-6_rk, 0.0_rk, 1.0e-6_rk]
        real(kind = rk), parameter :: A2(3) = [-1.30_rk, -1.33_rk, -1.36_rk]
        integer(kind = ik), parameter :: WANT(3) = [TQS_ERROR_NOT_STAR_CONVEX, &
                TQS_ERROR_NOT_STAR_CONVEX, TQS_ERROR_SEPARATED]
        real(kind = rk) :: radii(NT)
        type(result_t) :: res
        integer(kind = ik) :: j, status
        character(len = 40) :: tag
        do j = 1_ik, 3_ik
            write(tag, '(A,ES9.1)') 'K7 locus e = ', ES(j)
            p = pv(2.0_rk * sqrt(2.0_rk) + ES(j), -1.0_rk, 1.0_rk, 0.0_rk, 0.0_rk)
            call compute_radius_grid_standalone_s(p, thetas, radii, status)
            call assert_int_eq(status, WANT(j), trim(tag) // ': R-path status')
        end do
        do j = 1_ik, 3_ik
            write(tag, '(A,F6.2)') 'K7 pear alpha2 = ', A2(j)
            p = pv(1.0_rk, -1.0_rk, 1.0_rk, A2(j), 0.0_rk)
            call run_one_shot(p, thetas, 65_ik, res)
            call assert_true(all(res%status(1:N_OUT) == 0_ik), trim(tag) // ': every output valid')
        end do
    end subroutine test_k7_a3_through_zero

    subroutine test_k8_sphere_has_no_neck()
        real(kind = rk) :: zn, rn
        logical :: found
        integer(kind = ik) :: status
        call compute_neck_standalone_s(SPHERE_B, zn, rn, found, status)
        call assert_true(status == 0_ik .and. .not. found, 'K8 the sphere [0, -1, 1, 0, 0] has no neck')
    end subroutine test_k8_sphere_has_no_neck

    subroutine test_k9_large_coordinates()
        real(kind = rk) :: p(5)
        type(geometry_t) :: g
        integer(kind = ik) :: status
        real(kind = rk) :: i0, i1, i2, i4
        p = pv(1.0e8_rk, -1.0e-8_rk, 1.0_rk, 0.0_rk, 0.0_rk)
        call resolve_s(p, g, status)
        call assert_int_eq(status, 0_ik, 'K9 [1e8, -1e-8, 1, 0, 0] resolves')
        call assert_true(g%is_separated, 'K9 separated')
        call body_integrals_s(g, i0, i1, i2, i4)
        write(*, '(A,2ES12.3)') 'K9 volume rel. error, |z_cm|: ', abs(PI_C * i0 / (4.0_rk * PI_C / 3.0_rk) - 1.0_rk), &
                abs(i1 / i0)
        call assert_true(abs(PI_C * i0 / (4.0_rk * PI_C / 3.0_rk) - 1.0_rk) <= 1.0e-6_rk, 'K9 volume to 1e-6')
        call assert_true(abs(i1 / i0) <= 1.0e-6_rk, 'K9 centre of mass to 1e-6')
    end subroutine test_k9_large_coordinates

    subroutine test_k10_star_margin_edge()
        real(kind = rk) :: p(5)
        type(geometry_t) :: g
        real(kind = rk) :: lo, hi, mid, radii(NT)
        integer(kind = ik) :: j, status
        lo = 1.6_rk
        hi = 3.0_rk
        p = pv(lo, -0.6_rk, 1.0_rk, 0.9_rk, 0.0_rk)
        call resolve_s(p, g, status)
        call assert_true(status == 0_ik .and. star_margin_f(g) > STAR_CONVEXITY_MARGIN, 'K10 margin passes at 1.6')
        p = pv(hi, -0.6_rk, 1.0_rk, 0.9_rk, 0.0_rk)
        call resolve_s(p, g, status)
        call assert_true(status == 0_ik .and. star_margin_f(g) < STAR_CONVEXITY_MARGIN, 'K10 margin fails at 3.0')
        do j = 1_ik, 20_ik
            mid = 0.5_rk * (lo + hi)
            p = pv(mid, -0.6_rk, 1.0_rk, 0.9_rk, 0.0_rk)
            call resolve_s(p, g, status)
            if (star_margin_f(g) >= STAR_CONVEXITY_MARGIN) then
                lo = mid
            else
                hi = mid
            end if
        end do
        mid = 0.5_rk * (lo + hi)
        write(*, '(A,F12.8)') 'K10 star-margin edge on [s1, -0.6, 1, 0.9, 0] at s1 = ', mid
        p = pv(mid - 1.0e-3_rk, -0.6_rk, 1.0_rk, 0.9_rk, 0.0_rk)
        call compute_radius_grid_standalone_s(p, &
                thetas, radii, status)
        call assert_int_eq(status, 0_ik, 'K10 radius_grid just below the edge -> 0')
        p = pv(mid + 1.0e-3_rk, -0.6_rk, 1.0_rk, 0.9_rk, 0.0_rk)
        call compute_radius_grid_standalone_s(p, &
                thetas, radii, status)
        call assert_int_eq(status, TQS_ERROR_NOT_STAR_CONVEX, 'K10 radius_grid just above the edge -> 102')
    end subroutine test_k10_star_margin_edge

    subroutine test_k11_both_pass_count()
        real(kind = rk) :: p(5)
        type(geometry_t) :: g
        integer(kind = ik) :: i1, i2, i3, i4, i5, status, n_valid, n_both, n_points, n_ok
        n_both = 0_ik
        n_points = 0_ik
        n_ok = 0_ik
        do i1 = 0_ik, 14_ik
            do i2 = 0_ik, 12_ik
                do i3 = 0_ik, 6_ik
                    do i4 = 0_ik, 8_ik
                        do i5 = 0_ik, 8_ik
                            p = pv(0.25_rk * real(i1, rk), -2.0_rk + 0.25_rk * real(i2, rk), &
                                    0.5_rk + 0.25_rk * real(i3, rk), -0.8_rk + 0.2_rk * real(i4, rk), &
                                    -0.8_rk + 0.2_rk * real(i5, rk))
                            call resolve_unit_s(p, g, status, n_valid)
                            n_points = n_points + 1_ik
                            if (status == 0_ik) n_ok = n_ok + 1_ik
                            if (n_valid == 2_ik) n_both = n_both + 1_ik
                        end do
                    end do
                end do
            end do
        end do
        write(*, '(A,I0,A,I0,A,I0,A)') 'K11 both-pass roots: ', n_both, ' of ', n_points, &
                ' points (valid ', n_ok, ')'
        call assert_int_eq(n_points, 110565_ik, 'K11 box size')
    end subroutine test_k11_both_pass_count

end program tqs_param_breakdown_test
