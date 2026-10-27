!> Workers kernels: factored segment integrals against Gauss-Legendre
!! quadrature (polynomial integrands of degree <= 5, so 64 nodes are exact to
!! roundoff), the large-offset case the difference-of-primitives form loses,
!! body sums over a geometry record, and the stable quadratic.
program tqs_param_workers_test
    use precision_utilities_mod, only: ik, rk
    use mathematical_utilities_mod, only: compute_gauss_legendre_quadrature_s
    use tqs_workers_mod, only: geometry_t, segment_integrals_s, body_integrals_s, &
            solve_quadratic_s
    use test_utils_mod, only: assert_true, assert_int_eq, assert_close, &
            assert_abs_close, test_summary
    implicit none

    integer(kind = ik), parameter :: N_GL = 64_ik
    real(kind = rk) :: gl_x(N_GL), gl_w(N_GL)

    call compute_gauss_legendre_quadrature_s(N_GL, gl_x, gl_w)

    call test_w1_against_quadrature()
    call test_w2_large_offset()
    call test_w3_body_sum()
    call test_w4_quadratic()

    call test_summary()

contains

    !> The four integrands of segment_integrals_s by 64-node GL on [z_lo, z_hi].
    subroutine gl_segment(a_sq, k, l, z_lo, z_hi, i0, i1, i2, i4)
        real(kind = rk), intent(in) :: a_sq, k, l, z_lo, z_hi
        real(kind = rk), intent(out) :: i0, i1, i2, i4
        real(kind = rk) :: half, mid, z, rho_sq, wz
        integer(kind = ik) :: j
        half = 0.5_rk * (z_hi - z_lo)
        mid = 0.5_rk * (z_hi + z_lo)
        i0 = 0.0_rk
        i1 = 0.0_rk
        i2 = 0.0_rk
        i4 = 0.0_rk
        do j = 1_ik, N_GL
            z = mid + half * gl_x(j)
            wz = half * gl_w(j)
            rho_sq = a_sq - k * (z - l)**2
            i0 = i0 + wz * rho_sq
            i1 = i1 + wz * z * rho_sq
            i2 = i2 + wz * z * z * rho_sq
            i4 = i4 + wz * rho_sq * rho_sq
        end do
    end subroutine gl_segment

    subroutine test_w1_against_quadrature()
        real(kind = rk), parameter :: SEG(5, 3) = reshape([ &
                1.0_rk, 1.0_rk, 0.0_rk, -1.0_rk, 0.4_rk, &
                0.8_rk, -0.6_rk, 0.3_rk, -0.2_rk, 0.7_rk, &
                1.25_rk, 2.0_rk, -0.25_rk, -0.86_rk, -0.75_rk], [5, 3])
        real(kind = rk) :: i0, i1, i2, i4, g0, g1, g2, g4
        integer(kind = ik) :: s
        character(len = 16) :: tag
        do s = 1_ik, 3_ik
            write(tag, '(A,I0)') 'W1 seg ', s
            call segment_integrals_s(SEG(1, s), SEG(2, s), SEG(3, s), SEG(4, s), &
                    SEG(5, s), i0, i1, i2, i4)
            call gl_segment(SEG(1, s), SEG(2, s), SEG(3, s), SEG(4, s), SEG(5, s), &
                    g0, g1, g2, g4)
            call assert_abs_close(i0, g0, 1.0e-13_rk, trim(tag) // ' i0')
            call assert_abs_close(i1, g1, 1.0e-13_rk, trim(tag) // ' i1')
            call assert_abs_close(i2, g2, 1.0e-13_rk, trim(tag) // ' i2')
            call assert_abs_close(i4, g4, 1.0e-13_rk, trim(tag) // ' i4')
        end do
    end subroutine test_w1_against_quadrature

    subroutine test_w2_large_offset()
        real(kind = rk), parameter :: L_BIG = 5.0e7_rk
        real(kind = rk) :: i0, i1, i2, i4
        call segment_integrals_s(1.0_rk, 1.0_rk, L_BIG, L_BIG - 1.0_rk, L_BIG + 1.0_rk, &
                i0, i1, i2, i4)
        call assert_close(i0, 4.0_rk / 3.0_rk, 1.0e-12_rk, 'W2 i0 = 4/3')
        call assert_close(i1, L_BIG * 4.0_rk / 3.0_rk, 1.0e-12_rk, 'W2 i1 = l i0')
        call assert_close(i2, 4.0_rk / 15.0_rk + 2.5e15_rk * 4.0_rk / 3.0_rk, 1.0e-12_rk, &
                'W2 i2')
        call assert_close(i4, 16.0_rk / 15.0_rk, 1.0e-12_rk, 'W2 i4 = 16/15')
    end subroutine test_w2_large_offset

    subroutine test_w3_body_sum()
        type(geometry_t) :: geom
        real(kind = rk) :: b0, b1, b2, b4, s0, s1, s2, s4, p0, p1, p2, p4
        real(kind = rk), parameter :: W = 0.1_rk

        geom%a_sq = [1.0_rk, 1.0_rk, 0.47_rk]
        geom%k = [1.0_rk, 1.0_rk, -0.6_rk]
        geom%l = [-0.8_rk, 0.8_rk, 0.0_rk]
        geom%z_join = [-0.39_rk, 0.39_rk]
        geom%z_south = -1.8_rk
        geom%z_north = 1.8_rk

        call body_integrals_s(geom, b0, b1, b2, b4)
        call segment_integrals_s(geom%a_sq(1), geom%k(1), geom%l(1), geom%z_south, &
                geom%z_join(1), s0, s1, s2, s4)
        call segment_integrals_s(geom%a_sq(3), geom%k(3), geom%l(3), geom%z_join(1), &
                geom%z_join(2), p0, p1, p2, p4)
        s0 = s0 + p0; s1 = s1 + p1; s2 = s2 + p2; s4 = s4 + p4
        call segment_integrals_s(geom%a_sq(2), geom%k(2), geom%l(2), geom%z_join(2), &
                geom%z_north, p0, p1, p2, p4)
        s0 = s0 + p0; s1 = s1 + p1; s2 = s2 + p2; s4 = s4 + p4
        call assert_abs_close(b0, s0, 1.0e-15_rk, 'W3 connected i0 = segment sum')
        call assert_abs_close(b1, s1, 1.0e-15_rk, 'W3 connected i1 = segment sum')
        call assert_abs_close(b2, s2, 1.0e-15_rk, 'W3 connected i2 = segment sum')
        call assert_abs_close(b4, s4, 1.0e-15_rk, 'W3 connected i4 = segment sum')

        ! Separated: the middle contributes only outside the void (l3 - w, l3 + w).
        geom%is_separated = .true.
        geom%w = W
        call body_integrals_s(geom, b0, b1, b2, b4)
        call segment_integrals_s(geom%a_sq(1), geom%k(1), geom%l(1), geom%z_south, &
                geom%z_join(1), s0, s1, s2, s4)
        call segment_integrals_s(geom%a_sq(3), geom%k(3), geom%l(3), geom%z_join(1), &
                geom%l(3) - W, p0, p1, p2, p4)
        s0 = s0 + p0; s1 = s1 + p1; s2 = s2 + p2; s4 = s4 + p4
        call segment_integrals_s(geom%a_sq(3), geom%k(3), geom%l(3), geom%l(3) + W, &
                geom%z_join(2), p0, p1, p2, p4)
        s0 = s0 + p0; s1 = s1 + p1; s2 = s2 + p2; s4 = s4 + p4
        call segment_integrals_s(geom%a_sq(2), geom%k(2), geom%l(2), geom%z_join(2), &
                geom%z_north, p0, p1, p2, p4)
        s0 = s0 + p0; s1 = s1 + p1; s2 = s2 + p2; s4 = s4 + p4
        call assert_abs_close(b0, s0, 1.0e-15_rk, 'W3 separated i0 = segment sum')
        call assert_abs_close(b1, s1, 1.0e-15_rk, 'W3 separated i1 = segment sum')
        call assert_abs_close(b2, s2, 1.0e-15_rk, 'W3 separated i2 = segment sum')
        call assert_abs_close(b4, s4, 1.0e-15_rk, 'W3 separated i4 = segment sum')
    end subroutine test_w3_body_sum

    subroutine test_w4_quadratic()
        integer(kind = ik) :: n
        real(kind = rk) :: r1, r2

        call solve_quadratic_s(1.0_rk, -3.0_rk, 2.0_rk, n, r1, r2)
        call assert_int_eq(n, 2_ik, 'W4 (1,-3,2) two roots')
        call assert_abs_close(min(r1, r2), 1.0_rk, 1.0e-15_rk, 'W4 (1,-3,2) root 1')
        call assert_abs_close(max(r1, r2), 2.0_rk, 1.0e-15_rk, 'W4 (1,-3,2) root 2')

        call solve_quadratic_s(0.0_rk, 2.0_rk, -4.0_rk, n, r1, r2)
        call assert_int_eq(n, 1_ik, 'W4 (0,2,-4) linear: one root')
        call assert_abs_close(r1, 2.0_rk, 1.0e-15_rk, 'W4 (0,2,-4) root 2')

        call solve_quadratic_s(0.0_rk, 0.0_rk, 1.0_rk, n, r1, r2)
        call assert_int_eq(n, 0_ik, 'W4 (0,0,1) no root')

        call solve_quadratic_s(0.0_rk, 0.0_rk, 0.0_rk, n, r1, r2)
        call assert_int_eq(n, -1_ik, 'W4 (0,0,0) identically satisfied')

        call solve_quadratic_s(1.0_rk, 2.0_rk, 1.0_rk, n, r1, r2)
        call assert_int_eq(n, 1_ik, 'W4 (1,2,1) double root')
        call assert_abs_close(r1, -1.0_rk, 1.0e-15_rk, 'W4 (1,2,1) root -1')

        call solve_quadratic_s(1.0_rk, 2.0_rk, 1.0_rk + 1.0e-15_rk, n, r1, r2)
        call assert_int_eq(n, 1_ik, 'W4 (1,2,1+1e-15) clamped to a double root')
        call assert_abs_close(r1, -1.0_rk, 1.0e-15_rk, 'W4 clamped root -1')

        call solve_quadratic_s(1.0_rk, 2.0_rk, 1.01_rk, n, r1, r2)
        call assert_int_eq(n, 0_ik, 'W4 (1,2,1.01) no real root')

        call solve_quadratic_s(1.0e-8_rk, 1.0_rk, 1.0_rk, n, r1, r2)
        call assert_int_eq(n, 2_ik, 'W4 (1e-8,1,1) two roots')
        call assert_true(abs(1.0e-8_rk * r1**2 + r1 + 1.0_rk) <= 1.0e-12_rk * (1.0_rk + abs(r1)), &
                'W4 (1e-8,1,1) root 1 residual')
        call assert_true(abs(1.0e-8_rk * r2**2 + r2 + 1.0_rk) <= 1.0e-12_rk * (1.0_rk + abs(r2)), &
                'W4 (1e-8,1,1) root 2 residual')
    end subroutine test_w4_quadratic

end program tqs_param_workers_test
