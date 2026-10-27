!> Coordinates -> geometry for the three-quadratic-surface family, closed form.
!!
!! Nix's five coordinates [sigma1, sigma2, sigma3, alpha2, alpha3] fix the
!! two end spheroids directly (step 1). The two smooth-join conditions then
!! reduce to a quadratic in the waist position l3 (step 2), whose roots are
!! screened by the validity of the joins (step 3). The middle body is
!! classified by the sign of sigma2 and a3^2 and by where its waist lies. The
!! shape is finally scaled to unit volume (R0 = 1) and shifted to its centre
!! of mass.
!!
!! Degenerate curvature (an end with the middle's curvature, |k_i - sigma2|
!! at or below TQS_DEGENERATE_TOL) has no join of its own: the only
!! consistent limit is the two surfaces coinciding, accepted when the two end
!! surfaces meet at the shared join to within TQS_COINCIDENCE_TOL in rho^2.
module tqs_resolve_mod

    use precision_utilities_mod, only: ik, rk
    use shape_core_mod, only: SHAPE_VALID
    use tqs_workers_mod, only: geometry_t, solve_quadratic_s, body_integrals_s, TQS_PINCH_TOL, &
            TQS_DEGENERATE_TOL, TQS_COINCIDENCE_TOL, TQS_MAX_WAIST_OFFSET, TQS_MID_SPHEROID, &
            TQS_MID_CYLINDER, TQS_MID_ONE_SHEET, TQS_MID_CONE, TQS_MID_TWO_SHEET, &
            TQS_ERROR_INVALID_COORDINATES, TQS_ERROR_NO_GEOMETRY

    implicit none

    private

    !> z1 above z2 by at most this (times max(1, |z1|, |z2|), u units) is a
    !! zero-length middle, not a crossing: where the two end spheroids are
    !! tangent the l3 quadratic has a double root, z1 = z2 exactly, and
    !! roundoff (amplified by 1/|k_i - sigma2|) returns either order.
    real(kind = rk), parameter :: JOIN_ORDER_TOL = 1.0e-10_rk

    public :: validate_coordinates_s
    public :: resolve_unit_s
    public :: resolve_s
    public :: alpha1_f
    public :: scission_sigma1_s

contains

    !> Code 100 unless k1 = sigma3 + alpha3/2 > 0, k2 = sigma3 - alpha3/2 > 0,
    !! |alpha2| < 2 and sigma1 >= 0.
    !!
    !! @param[in]  params  [sigma1, sigma2, sigma3, alpha2, alpha3]
    !! @param[out] status  SHAPE_VALID or TQS_ERROR_INVALID_COORDINATES
    pure subroutine validate_coordinates_s(params, status)

        real(kind = rk), intent(in) :: params(:)
        integer(kind = ik), intent(out) :: status

        real(kind = rk) :: k1, k2

        k1 = params(3) + 0.5_rk * params(5)
        k2 = params(3) - 0.5_rk * params(5)
        if (k1 > 0.0_rk .and. k2 > 0.0_rk .and. abs(params(4)) < 2.0_rk &
                .and. params(1) >= 0.0_rk) then
            status = SHAPE_VALID
        else
            status = TQS_ERROR_INVALID_COORDINATES
        end if

    end subroutine validate_coordinates_s

    !> Coordinates -> geometry in units of u (u = 1), unshifted.
    !!
    !! Fills a_sq, k, c, l, z_join, z_south = l1 - c1, z_north = l2 + c2, w,
    !! mid_kind and the three flags. On a nonzero status the geometry is the
    !! zero-filled default record.
    !!
    !! @param[in]  params         [sigma1, sigma2, sigma3, alpha2, alpha3]
    !! @param[out] geom           The unit-stage geometry
    !! @param[out] status         SHAPE_VALID, 100 or 101
    !! @param[out] n_valid_roots  Generic-path roots that passed the validity
    !!                            checks (0, 1 or 2); 0 on the degenerate paths
    pure subroutine resolve_unit_s(params, geom, status, n_valid_roots)

        real(kind = rk), intent(in) :: params(:)
        type(geometry_t), intent(out) :: geom
        integer(kind = ik), intent(out) :: status
        integer(kind = ik), intent(out), optional :: n_valid_roots

        real(kind = rk) :: s1, s2, s3, a2, a3
        real(kind = rk) :: g(2), coef_a, coef_b, coef_c, r(2), cand(2), l_mid
        real(kind = rk) :: l3, a3_sq, z1, z2, best_l3, best_a3_sq, best_z1, best_z2
        integer(kind = ik) :: n_roots, n_cand, n_valid, j, i_deg, j_oth
        logical :: deg(2), ok

        if (present(n_valid_roots)) n_valid_roots = 0_ik

        call validate_coordinates_s(params, status)
        if (status /= SHAPE_VALID) return

        s1 = params(1)
        s2 = params(2)
        s3 = params(3)
        a2 = params(4)
        a3 = params(5)

        ! Step 1: the end spheroids.
        geom%u = 1.0_rk
        geom%k(1) = s3 + 0.5_rk * a3
        geom%k(2) = s3 - 0.5_rk * a3
        geom%k(3) = s2
        geom%a_sq(1) = 1.0_rk + 0.5_rk * a2
        geom%a_sq(2) = 1.0_rk - 0.5_rk * a2
        geom%c(1) = sqrt(geom%a_sq(1) / geom%k(1))
        geom%c(2) = sqrt(geom%a_sq(2) / geom%k(2))
        geom%l(1) = -0.5_rk * s1
        geom%l(2) = 0.5_rk * s1
        l_mid = 0.5_rk * (geom%l(1) + geom%l(2))

        deg(1) = abs(geom%k(1) - s2) <= TQS_DEGENERATE_TOL
        deg(2) = abs(geom%k(2) - s2) <= TQS_DEGENERATE_TOL

        if (deg(1) .and. deg(2)) then
            ! Both ends share the middle's curvature: only the spheroid line
            ! sigma1 = alpha2 = 0 is consistent, decided exactly.
            if (abs(s1) > 0.0_rk .or. abs(a2) > 0.0_rk) then
                call fail_s(geom, status)
                return
            end if
            l3 = 0.0_rk
            a3_sq = 1.0_rk
            z1 = 0.0_rk
            z2 = 0.0_rk
            if (.not. joins_valid_f(geom, z1, z2)) then
                call fail_s(geom, status)
                return
            end if

        else if (deg(1) .or. deg(2)) then
            ! One degenerate end i: its surface must coincide with the middle.
            ! l3 = l_i, a3^2 from the other join j; the coincident join is
            ! arbitrary on the shared surface and is set equal to join j, so
            ! the shape is end i meeting end j there. Coincidence is judged by
            ! how well they meet: the rho^2 mismatch of the two end surfaces at
            ! that join, which a close a3^2 alone does not bound.
            if (deg(1)) then
                i_deg = 1_ik
                j_oth = 2_ik
            else
                i_deg = 2_ik
                j_oth = 1_ik
            end if
            l3 = geom%l(i_deg)
            a3_sq = geom%a_sq(j_oth) + join_g_f(geom%k(j_oth), s2) * (l3 - geom%l(j_oth))**2
            z2 = join_z_f(geom%k(j_oth), geom%l(j_oth), s2, l3)
            z1 = z2
            if (abs(geom%a_sq(i_deg) - geom%k(i_deg) * (z2 - geom%l(i_deg))**2 &
                    - geom%a_sq(j_oth) + geom%k(j_oth) * (z2 - geom%l(j_oth))**2) &
                    > TQS_COINCIDENCE_TOL * max(1.0_rk, geom%a_sq(i_deg))) then
                call fail_s(geom, status)
                return
            end if
            if (.not. joins_valid_f(geom, z1, z2)) then
                call fail_s(geom, status)
                return
            end if

        else
            ! Generic path: equate the two joins' a3^2, a quadratic in l3.
            g(1) = join_g_f(geom%k(1), s2)
            g(2) = join_g_f(geom%k(2), s2)
            coef_a = g(1) - g(2)
            coef_b = -2.0_rk * (g(1) * geom%l(1) - g(2) * geom%l(2))
            coef_c = g(1) * geom%l(1)**2 - g(2) * geom%l(2)**2 + geom%a_sq(1) - geom%a_sq(2)
            call solve_quadratic_s(coef_a, coef_b, coef_c, n_roots, r(1), r(2))

            select case (n_roots)
            case (-1_ik)
                n_cand = 1_ik
                cand(1) = l_mid
            case (0_ik)
                call fail_s(geom, status)
                return
            case (1_ik)
                n_cand = 1_ik
                cand(1) = r(1)
            case default
                n_cand = 2_ik
                cand = r
            end select

            n_valid = 0_ik
            best_l3 = 0.0_rk
            best_a3_sq = 0.0_rk
            best_z1 = 0.0_rk
            best_z2 = 0.0_rk
            do j = 1_ik, n_cand
                l3 = cand(j)
                z1 = join_z_f(geom%k(1), geom%l(1), s2, l3)
                z2 = join_z_f(geom%k(2), geom%l(2), s2, l3)
                ok = joins_valid_f(geom, z1, z2) .and. abs(l3 - l_mid) <= TQS_MAX_WAIST_OFFSET
                if (.not. ok) cycle
                n_valid = n_valid + 1_ik
                ! Both pass: the root nearer the end-centre midpoint.
                if (n_valid > 1_ik) then
                    if (.not. (abs(l3 - l_mid) < abs(best_l3 - l_mid))) cycle
                end if
                best_l3 = l3
                best_a3_sq = geom%a_sq(1) + g(1) * (l3 - geom%l(1))**2
                best_z1 = z1
                best_z2 = z2
            end do
            if (present(n_valid_roots)) n_valid_roots = n_valid
            if (n_valid < 1_ik) then
                call fail_s(geom, status)
                return
            end if
            l3 = best_l3
            a3_sq = best_a3_sq
            z1 = best_z1
            z2 = best_z2
            ! A tangent double root: one join, the midpoint of the two roundings.
            if (z1 > z2) then
                z1 = 0.5_rk * (z1 + z2)
                z2 = z1
            end if
        end if

        geom%l(3) = l3
        geom%a_sq(3) = a3_sq
        geom%z_join(1) = z1
        geom%z_join(2) = z2
        geom%z_south = geom%l(1) - geom%c(1)
        geom%z_north = geom%l(2) + geom%c(2)

        call classify_s(geom)

    end subroutine resolve_unit_s

    !> Coordinates -> geometry in R0 = 1 units with the centre of mass at the origin.
    !!
    !! resolve_unit_s, then the volume scale u = (4 pi/3 / V1)^(1/3) from the
    !! unit-stage volume V1 = pi i0: lengths (c, l, z_join, z_south, z_north, w)
    !! times u, a_sq times u^2, k unchanged. The centre of mass z_cm = i1/i0,
    !! computed at u = 1 and scaled by u, is then subtracted from l, z_join,
    !! z_south and z_north. geom%u = u.
    !!
    !! @param[in]  params         [sigma1, sigma2, sigma3, alpha2, alpha3]
    !! @param[out] geom           The resolved geometry; zero-filled on failure
    !! @param[out] status         SHAPE_VALID, 100 or 101
    !! @param[out] n_valid_roots  As resolve_unit_s
    pure subroutine resolve_s(params, geom, status, n_valid_roots)

        real(kind = rk), intent(in) :: params(:)
        type(geometry_t), intent(out) :: geom
        integer(kind = ik), intent(out) :: status
        integer(kind = ik), intent(out), optional :: n_valid_roots

        real(kind = rk) :: i0, i1, i2, i4, u, shift

        call resolve_unit_s(params, geom, status, n_valid_roots)
        if (status /= SHAPE_VALID) return

        call body_integrals_s(geom, i0, i1, i2, i4)
        ! pi i0 = 4 pi/3 u^-3, so u^3 = 4/(3 i0).
        u = (4.0_rk / (3.0_rk * i0))**(1.0_rk / 3.0_rk)
        shift = u * (i1 / i0)

        geom%u = u
        geom%a_sq = geom%a_sq * (u * u)
        geom%c = geom%c * u
        geom%w = geom%w * u
        geom%l = geom%l * u - shift
        geom%z_join = geom%z_join * u - shift
        geom%z_south = geom%z_south * u - shift
        geom%z_north = geom%z_north * u - shift

    end subroutine resolve_s

    !> Nix's alpha1 = (l1 + l2)/(2u): where the centre of mass put the end-centre midpoint.
    pure function alpha1_f(geom) result(alpha1)
        type(geometry_t), intent(in) :: geom
        real(kind = rk) :: alpha1
        alpha1 = 0.5_rk * (geom%l(1) + geom%l(2)) / geom%u
    end function alpha1_f

    !> The scission map of (sigma2, sigma3, alpha2, alpha3), in units of u.
    !!
    !! sigma1_t = c1 + c2 with c_i = sqrt((2 +- alpha2)/(2 sigma3 +- alpha3)):
    !! Moller's eq. 16, the two end spheroids touching pole to pole. For
    !! sigma2 < 0 the a3^2 = 0 locus of the resolve, |l3 - l_i| =
    !! c_i (1 + k_i/|sigma2|)^(1/2), gives sigma1_star = c1 (1 + k1/|s2|)^(1/2)
    !! + c2 (1 + k2/|s2|)^(1/2); a spheroid or cylinder middle never pinches,
    !! so for sigma2 >= 0 sigma1_star = 0 and defined = .false. sigma1 does not
    !! enter; it is validated like every other output.
    !!
    !! @param[in]  params       [sigma1, sigma2, sigma3, alpha2, alpha3]
    !! @param[out] sigma1_t     Touching-spheroids separation
    !! @param[out] sigma1_star  Scission separation at this sigma2; 0 unless defined
    !! @param[out] defined      .true. iff sigma2 < 0 and status is 0
    !! @param[out] status       SHAPE_VALID or 100
    pure subroutine scission_sigma1_s(params, sigma1_t, sigma1_star, defined, status)

        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(out) :: sigma1_t, sigma1_star
        logical, intent(out) :: defined
        integer(kind = ik), intent(out) :: status

        real(kind = rk) :: s2, k1, k2, c1, c2

        sigma1_t = 0.0_rk
        sigma1_star = 0.0_rk
        defined = .false.

        call validate_coordinates_s(params, status)
        if (status /= SHAPE_VALID) return

        s2 = params(2)
        k1 = params(3) + 0.5_rk * params(5)
        k2 = params(3) - 0.5_rk * params(5)
        c1 = sqrt((2.0_rk + params(4)) / (2.0_rk * params(3) + params(5)))
        c2 = sqrt((2.0_rk - params(4)) / (2.0_rk * params(3) - params(5)))
        sigma1_t = c1 + c2

        if (s2 < 0.0_rk) then
            sigma1_star = c1 * sqrt(1.0_rk + k1 / abs(s2)) + c2 * sqrt(1.0_rk + k2 / abs(s2))
            defined = .true.
        end if

    end subroutine scission_sigma1_s

    !> g_i = k_i sigma2 / (k_i - sigma2): a3^2 = a_i^2 + g_i (l3 - l_i)^2 at join i.
    pure function join_g_f(k_end, s2) result(g)
        real(kind = rk), intent(in) :: k_end, s2
        real(kind = rk) :: g
        g = k_end * s2 / (k_end - s2)
    end function join_g_f

    !> Join position from the derivative condition k_i (z - l_i) = sigma2 (z - l3).
    pure function join_z_f(k_end, l_end, s2, l3) result(z)
        real(kind = rk), intent(in) :: k_end, l_end, s2, l3
        real(kind = rk) :: z
        z = l_end + s2 * (l_end - l3) / (k_end - s2)
    end function join_z_f

    !> z1 <= z2 (up to JOIN_ORDER_TOL) and each join strictly inside its end
    !! spheroid's z-extent.
    pure function joins_valid_f(geom, z1, z2) result(ok)
        type(geometry_t), intent(in) :: geom
        real(kind = rk), intent(in) :: z1, z2
        logical :: ok
        ok = z1 - z2 <= JOIN_ORDER_TOL * max(1.0_rk, abs(z1), abs(z2)) &
                .and. geom%l(1) - geom%c(1) < z1 .and. z1 < geom%l(1) + geom%c(1) &
                .and. geom%l(2) - geom%c(2) < z2 .and. z2 < geom%l(2) + geom%c(2)
    end function joins_valid_f

    !> Middle-body kind and the waist flags of a geometry with its joins set.
    pure subroutine classify_s(geom)

        type(geometry_t), intent(inout) :: geom

        real(kind = rk) :: s2, a3_sq

        s2 = geom%k(3)
        a3_sq = geom%a_sq(3)
        ! A waist exists only for sigma2 < 0 (a convex middle has an equator,
        ! not a waist). Strict: a zero-length middle or a waist on a join is
        ! never inside.
        geom%waist_inside = s2 < 0.0_rk .and. geom%z_join(1) < geom%l(3) &
                .and. geom%l(3) < geom%z_join(2)

        if (s2 > 0.0_rk) then
            geom%mid_kind = TQS_MID_SPHEROID
        else if (.not. (abs(s2) > 0.0_rk)) then
            geom%mid_kind = TQS_MID_CYLINDER
        else if (abs(a3_sq) <= TQS_PINCH_TOL) then
            ! Counted as zero, and stored as zero: the middle is the exact
            ! cone, so a pinch point is rho = 0 whatever sign roundoff gave.
            geom%mid_kind = TQS_MID_CONE
            geom%a_sq(3) = 0.0_rk
            geom%is_pinched = geom%waist_inside
        else if (a3_sq > 0.0_rk) then
            geom%mid_kind = TQS_MID_ONE_SHEET
        else
            geom%mid_kind = TQS_MID_TWO_SHEET
            geom%is_separated = geom%waist_inside
            if (geom%is_separated) geom%w = sqrt(a3_sq / s2)
        end if

    end subroutine classify_s

    !> The zero-filled record that every failure reports.
    pure subroutine fail_s(geom, status)
        type(geometry_t), intent(out) :: geom
        integer(kind = ik), intent(out) :: status
        status = TQS_ERROR_NO_GEOMETRY
        geom = geometry_t()
    end subroutine fail_s

end module tqs_resolve_mod
