!> Arithmetic on a resolved three-quadratic-surface (3QS) geometry.
!!
!! Owns the constants and enums of the library, the geometry record, and
!! every kernel that works on a resolved geometry: segment and body
!! integrals, the stable quadratic, the profile, the neck, ray intersection,
!! the star-convexity margin and the Moller coordinates. Everything is closed
!! form: no Newton iteration, no quadrature, no root search.
module tqs_workers_mod

    use precision_utilities_mod, only: ik, rk

    implicit none

    private

    !---------------------------------------------------------------------------
    ! Fixed arity and resolution floor
    !---------------------------------------------------------------------------
    !> params = [sigma1, sigma2, sigma3, alpha2, alpha3]; positions are permanent.
    integer(kind = ik), parameter, public :: TQS_N_PARAMS = 5_ik
    !> The two poles and one interior node: a sanity bound, the profile is closed form.
    integer(kind = ik), parameter, public :: TQS_N_POINTS_FLOOR = 3_ik

    !---------------------------------------------------------------------------
    ! Tolerances
    !---------------------------------------------------------------------------
    !> Minimum of rho - z drho/dz for the R(theta) outputs (representability, as fos-param).
    real(kind = rk), parameter, public :: STAR_CONVEXITY_MARGIN = 1.0e-2_rk
    !> |a3^2| in u^2 units at or below this counts as zero: roundoff at the
    !! scission locus must not flip the cone / two_sheet split.
    real(kind = rk), parameter, public :: TQS_PINCH_TOL = 1.0e-12_rk
    !> |k_i - sigma2| at or below this counts as equal (degenerate curvature).
    real(kind = rk), parameter, public :: TQS_DEGENERATE_TOL = 1.0e-4_rk
    !> A degenerate end coincides with the middle when the rho^2 mismatch of the
    !! two end surfaces at their shared join is at most tol * max(1, a_i^2).
    !! On a true (tangent) coincident family the mismatch is O((k_i - sigma2)^2),
    !! below 1e-7 inside TQS_DEGENERATE_TOL; surfaces that merely have close
    !! a_i^2 leave a visible step in rho.
    real(kind = rk), parameter, public :: TQS_COINCIDENCE_TOL = 1.0e-6_rk
    !> |l3 - (l1 + l2)/2| above this (u units) is code 101. Near sigma2 = 0 with
    !! alpha2 /= 0 the waist runs off like alpha2/sigma2 and a3^2 like 1/sigma2;
    !! the middle a3^2 - sigma2 (z - l3)^2 is then a difference of huge terms
    !! and the centre of mass loses ~eps * offset^2. At 1e3 that is ~1e-11 R0;
    !! the rejected band is |sigma2| below ~1e-4, around the sigma2 = 0 point
    !! that is already 101 when alpha2 /= 0.
    real(kind = rk), parameter, public :: TQS_MAX_WAIST_OFFSET = 1.0e3_rk
    !> A discriminant negative by less than this fraction of its scale is a
    !! double root: on the coincident family the exact discriminant is zero.
    real(kind = rk), parameter :: QUADRATIC_DISC_TOL = 1.0e-12_rk
    !> Slack on a ray hit's z against a join: the two segments meeting at a
    !! join each round their own hit, and a ray exactly at the join angle must
    !! not fall between them. Lengths are of order 1 (R0 units).
    real(kind = rk), parameter :: RAY_JOIN_SLACK = 1.0e-12_rk

    !---------------------------------------------------------------------------
    ! Middle-body kinds
    !---------------------------------------------------------------------------
    integer(kind = ik), parameter, public :: TQS_MID_NONE = 0_ik
    integer(kind = ik), parameter, public :: TQS_MID_SPHEROID = 1_ik
    integer(kind = ik), parameter, public :: TQS_MID_CYLINDER = 2_ik
    integer(kind = ik), parameter, public :: TQS_MID_ONE_SHEET = 3_ik
    integer(kind = ik), parameter, public :: TQS_MID_CONE = 4_ik
    integer(kind = ik), parameter, public :: TQS_MID_TWO_SHEET = 5_ik

    !---------------------------------------------------------------------------
    ! Library status codes (100+); the shared 0-5 live in shape_core_mod
    !---------------------------------------------------------------------------
    !! Owned here because the resolve module raises them and cannot depend on
    !! the front module; the front module re-exports them.
    integer(kind = ik), parameter, public :: TQS_ERROR_INVALID_COORDINATES = 100_ik
    integer(kind = ik), parameter, public :: TQS_ERROR_NO_GEOMETRY = 101_ik
    integer(kind = ik), parameter, public :: TQS_ERROR_NOT_STAR_CONVEX = 102_ik
    integer(kind = ik), parameter, public :: TQS_ERROR_SEPARATED = 103_ik

    !---------------------------------------------------------------------------
    ! Geometry record
    !---------------------------------------------------------------------------
    !> A resolved shape. Segment i is the surface rho^2 = a_sq(i) - k(i) (z - l(i))^2:
    !! 1 the south end spheroid on [z_south, z_join(1)], 3 the middle body on
    !! [z_join(1), z_join(2)], 2 the north end spheroid on [z_join(2), z_north].
    !! c(i) = sqrt(a_sq(i)/k(i)) are the end semi-axes along z. Lengths are in
    !! units of u until the resolve rescales them; u is Nix's scale.
    type, public :: geometry_t
        real(kind = rk) :: u = 0.0_rk
        real(kind = rk) :: a_sq(3) = 0.0_rk
        real(kind = rk) :: k(3) = 0.0_rk
        real(kind = rk) :: c(2) = 0.0_rk
        real(kind = rk) :: l(3) = 0.0_rk
        real(kind = rk) :: z_join(2) = 0.0_rk
        real(kind = rk) :: z_south = 0.0_rk
        real(kind = rk) :: z_north = 0.0_rk
        real(kind = rk) :: w = 0.0_rk            !! Half-width of the void, 0 unless is_separated
        integer(kind = ik) :: mid_kind = TQS_MID_NONE
        logical :: waist_inside = .false.
        logical :: is_separated = .false.
        logical :: is_pinched = .false.
    end type geometry_t

    public :: segment_integrals_s
    public :: body_integrals_s
    public :: solve_quadratic_s
    public :: segment_rho_s
    public :: rho_at_z_s
    public :: profile_s
    public :: neck_s
    public :: poles_s
    public :: moller_coordinates_s
    public :: ray_radius_s
    public :: star_margin_f
    public :: theta_trig_s

contains

    !===========================================================================
    ! INTEGRALS
    !===========================================================================

    !> Moments of one surface rho^2 = a_sq - k (z - l)^2 over [z_lo, z_hi].
    !!
    !! Every integral is written as (z_hi - z_lo) * P(t_lo, t_hi) with
    !! t = z - l, never as a difference of primitives F(t_hi) - F(t_lo): the
    !! difference cancels catastrophically when |t| is large (sigma1 = 1e8
    !! loses every digit), the factored form does not. The moments in z follow
    !! from the moments in t by the binomial shift, which only multiplies.
    !!
    !! @param[in]  a_sq, k, l  The surface
    !! @param[in]  z_lo, z_hi  The interval; z_hi < z_lo gives the negated integrals
    !! @param[out] i0          int rho^2 dz
    !! @param[out] i1          int z rho^2 dz
    !! @param[out] i2          int z^2 rho^2 dz
    !! @param[out] i4          int rho^4 dz
    pure subroutine segment_integrals_s(a_sq, k, l, z_lo, z_hi, i0, i1, i2, i4)

        real(kind = rk), intent(in) :: a_sq, k, l, z_lo, z_hi
        real(kind = rk), intent(out) :: i0, i1, i2, i4

        real(kind = rk) :: d, ta, tb, s2, s4, j0, j1, j2

        d = z_hi - z_lo
        ta = z_lo - l
        tb = z_hi - l
        s2 = tb * tb + tb * ta + ta * ta
        s4 = tb**4 + tb**3 * ta + tb * tb * ta * ta + tb * ta**3 + ta**4

        j0 = d * (a_sq - k * s2 / 3.0_rk)
        j1 = d * (0.5_rk * a_sq * (tb + ta) - 0.25_rk * k * (tb + ta) * (tb * tb + ta * ta))
        j2 = d * (a_sq * s2 / 3.0_rk - k * s4 / 5.0_rk)

        i0 = j0
        i1 = j1 + l * j0
        i2 = j2 + 2.0_rk * l * j1 + l * l * j0
        i4 = d * (a_sq * a_sq - 2.0_rk * a_sq * k * s2 / 3.0_rk + k * k * s4 / 5.0_rk)

    end subroutine segment_integrals_s

    !> Moments of the whole body: end 1, the middle (minus the void when
    !! separated), end 2, summed in that order.
    !!
    !! @param[in]  geom            Resolved geometry, any units
    !! @param[out] i0, i1, i2, i4  As segment_integrals_s, over the body
    pure subroutine body_integrals_s(geom, i0, i1, i2, i4)

        type(geometry_t), intent(in) :: geom
        real(kind = rk), intent(out) :: i0, i1, i2, i4

        real(kind = rk) :: p0, p1, p2, p4

        call segment_integrals_s(geom%a_sq(1), geom%k(1), geom%l(1), geom%z_south, &
                geom%z_join(1), i0, i1, i2, i4)

        if (geom%is_separated) then
            call segment_integrals_s(geom%a_sq(3), geom%k(3), geom%l(3), geom%z_join(1), &
                    geom%l(3) - geom%w, p0, p1, p2, p4)
            call accumulate_s(p0, p1, p2, p4, i0, i1, i2, i4)
            call segment_integrals_s(geom%a_sq(3), geom%k(3), geom%l(3), &
                    geom%l(3) + geom%w, geom%z_join(2), p0, p1, p2, p4)
        else
            call segment_integrals_s(geom%a_sq(3), geom%k(3), geom%l(3), geom%z_join(1), &
                    geom%z_join(2), p0, p1, p2, p4)
        end if
        call accumulate_s(p0, p1, p2, p4, i0, i1, i2, i4)

        call segment_integrals_s(geom%a_sq(2), geom%k(2), geom%l(2), geom%z_join(2), &
                geom%z_north, p0, p1, p2, p4)
        call accumulate_s(p0, p1, p2, p4, i0, i1, i2, i4)

    end subroutine body_integrals_s

    pure subroutine accumulate_s(p0, p1, p2, p4, i0, i1, i2, i4)
        real(kind = rk), intent(in) :: p0, p1, p2, p4
        real(kind = rk), intent(inout) :: i0, i1, i2, i4
        i0 = i0 + p0
        i1 = i1 + p1
        i2 = i2 + p2
        i4 = i4 + p4
    end subroutine accumulate_s

    !===========================================================================
    ! QUADRATIC
    !===========================================================================

    !> Real roots of a x^2 + b x + c = 0 in the stable form.
    !!
    !! q = -(b + sign(b) sqrt(disc))/2, roots q/a and c/q, so neither root is
    !! a difference of nearly equal numbers. A discriminant negative by less
    !! than QUADRATIC_DISC_TOL * max(b^2, |4 a c|) is clamped to zero: on the
    !! coincident family the exact discriminant is zero for every parameter
    !! value and roundoff gives either sign.
    !!
    !! @param[in]  a, b, c  Coefficients
    !! @param[out] n_roots  2; 1 (double root, or linear when a = 0); 0 (no real
    !!                      root, or a = b = 0 with c /= 0); -1 (a = b = c = 0,
    !!                      identically satisfied)
    !! @param[out] r1, r2   The roots; r2 = r1 when n_roots = 1, both 0 when n_roots <= 0
    pure subroutine solve_quadratic_s(a, b, c, n_roots, r1, r2)

        real(kind = rk), intent(in) :: a, b, c
        integer(kind = ik), intent(out) :: n_roots
        real(kind = rk), intent(out) :: r1, r2

        real(kind = rk) :: disc, scale, q

        r1 = 0.0_rk
        r2 = 0.0_rk

        if (.not. (abs(a) > 0.0_rk)) then
            if (abs(b) > 0.0_rk) then
                n_roots = 1_ik
                r1 = -c / b
                r2 = r1
            else if (abs(c) > 0.0_rk) then
                n_roots = 0_ik
            else
                n_roots = -1_ik
            end if
            return
        end if

        disc = b * b - 4.0_rk * a * c
        if (disc > 0.0_rk) then
            q = -0.5_rk * (b + sign(sqrt(disc), b))
            n_roots = 2_ik
            r1 = q / a
            r2 = c / q
            return
        end if

        scale = max(b * b, abs(4.0_rk * a * c))
        if (-disc <= QUADRATIC_DISC_TOL * scale) then
            n_roots = 1_ik
            r1 = -0.5_rk * b / a
            r2 = r1
        else
            n_roots = 0_ik
        end if

    end subroutine solve_quadratic_s

    !===========================================================================
    ! PROFILE
    !===========================================================================

    !> One surface: rho = sqrt(a_sq - k (z - l)^2), drho/dz = -k (z - l)/rho.
    !!
    !! (0, 0) where rho^2 <= 0: outside the surface's extent, inside a void,
    !! and at a pole or pinch point, where drho/dz diverges (the fos-param
    !! convention).
    pure subroutine segment_rho_s(a_sq, k, l, z, rho, drho_dz)

        real(kind = rk), intent(in) :: a_sq, k, l, z
        real(kind = rk), intent(out) :: rho, drho_dz

        real(kind = rk) :: rho_sq

        rho_sq = a_sq - k * (z - l)**2
        if (rho_sq > 0.0_rk) then
            rho = sqrt(rho_sq)
            drho_dz = -k * (z - l) / rho
        else
            rho = 0.0_rk
            drho_dz = 0.0_rk
        end if

    end subroutine segment_rho_s

    !> rho and drho/dz of the body at z: end 1 for z <= z_join(1), end 2 for
    !! z >= z_join(2), the middle between.
    pure subroutine rho_at_z_s(geom, z, rho, drho_dz)

        type(geometry_t), intent(in) :: geom
        real(kind = rk), intent(in) :: z
        real(kind = rk), intent(out) :: rho, drho_dz

        integer(kind = ik) :: seg

        if (z <= geom%z_join(1)) then
            seg = 1_ik
        else if (z >= geom%z_join(2)) then
            seg = 2_ik
        else
            seg = 3_ik
        end if
        call segment_rho_s(geom%a_sq(seg), geom%k(seg), geom%l(seg), z, rho, drho_dz)

    end subroutine rho_at_z_s

    !> The profile on size(z) uniform nodes over [z_south, z_north].
    !!
    !! Node j sits at z_south + (z_north - z_south) * (j - 1)/(n - 1), the
    !! fraction formed first so a symmetric shape puts its middle node at
    !! exactly 0; z(1) and z(n) are the poles exactly, with (0, 0) there.
    !!
    !! @param[in]  geom     Resolved geometry
    !! @param[out] z        Nodes; size >= 2 (the caller checks the floor)
    !! @param[out] rho      rho(z), size(z) long
    !! @param[out] drho_dz  drho/dz, size(z) long
    pure subroutine profile_s(geom, z, rho, drho_dz)

        type(geometry_t), intent(in) :: geom
        real(kind = rk), intent(out) :: z(:), rho(:), drho_dz(:)

        integer(kind = ik) :: n, j
        real(kind = rk) :: span

        n = size(z, kind = ik)
        span = geom%z_north - geom%z_south
        z(1) = geom%z_south
        rho(1) = 0.0_rk
        drho_dz(1) = 0.0_rk
        do j = 2_ik, n - 1_ik
            z(j) = geom%z_south + span * (real(j - 1_ik, rk) / real(n - 1_ik, rk))
            call rho_at_z_s(geom, z(j), rho(j), drho_dz(j))
        end do
        z(n) = geom%z_north
        rho(n) = 0.0_rk
        drho_dz(n) = 0.0_rk

    end subroutine profile_s

    !===========================================================================
    ! NECK AND POLES
    !===========================================================================

    !> The waist of a middle body with sigma2 < 0 whose waist lies strictly
    !! inside it: z_neck = l3; rho_neck = a3 for a one-sheet middle, 0 for a
    !! pinched or separated one (the split plane). Otherwise found = .false.
    !! and both outputs are 0.
    pure subroutine neck_s(geom, z_neck, rho_neck, found)

        type(geometry_t), intent(in) :: geom
        real(kind = rk), intent(out) :: z_neck, rho_neck
        logical, intent(out) :: found

        z_neck = 0.0_rk
        rho_neck = 0.0_rk
        found = geom%k(3) < 0.0_rk .and. geom%waist_inside
        if (.not. found) return

        z_neck = geom%l(3)
        if (geom%mid_kind == TQS_MID_ONE_SHEET) rho_neck = sqrt(geom%a_sq(3))

    end subroutine neck_s

    !> Pole radii from the centre of mass: r_north = z_north, r_south = -z_south.
    pure subroutine poles_s(geom, r_north, r_south)
        type(geometry_t), intent(in) :: geom
        real(kind = rk), intent(out) :: r_north, r_south
        r_north = geom%z_north
        r_south = -geom%z_south
    end subroutine poles_s

    !===========================================================================
    ! R(THETA) AND THE STAR GATE
    !===========================================================================

    !> cos and sin of a theta set. The cache's primary set and a caller's
    !! at-thetas set both go through this one kernel, so the same thetas give
    !! the same bits on either path.
    pure subroutine theta_trig_s(thetas, cos_t, sin_t)
        real(kind = rk), intent(in) :: thetas(:)
        real(kind = rk), intent(out) :: cos_t(:), sin_t(:)
        integer(kind = ik) :: i
        do i = 1_ik, size(thetas, kind = ik)
            cos_t(i) = cos(thetas(i))
            sin_t(i) = sin(thetas(i))
        end do
    end subroutine theta_trig_s

    !> R(theta) from the origin by closed-form ray intersection.
    !!
    !! The ray z = r cos, rho = r sin on segment i satisfies
    !!   F(r, theta) = r^2 (sin^2 + k cos^2) - 2 k l cos r + (k l^2 - a_sq) = 0,
    !! solved in the stable form; the leading coefficient vanishes where the
    !! ray is parallel to a hyperboloid's asymptotic cone and the equation is
    !! then linear (solve_quadratic_s handles both). A root counts when r > 0
    !! and z = r cos lies on the segment's side of its joins (an end spheroid's
    !! outer bound, its pole, holds automatically on the surface). The first
    !! survivor in the order end 1, middle, end 2 is R; for a star-convex shape
    !! it is the only one. dR/dtheta = -F_theta/F_r with
    !!   F_r = 2 r (sin^2 + k cos^2) - 2 k l cos,
    !!   F_theta = 2 r sin [r cos (1 - k) + k l].
    !!
    !! @param[in]  geom                   Resolved, connected geometry
    !! @param[in]  cos_theta, sin_theta   The ray direction
    !! @param[out] r, dr_dtheta           R and dR/dtheta; 0 when not found
    !! @param[out] found                  .false. if no segment is hit
    pure subroutine ray_radius_s(geom, cos_theta, sin_theta, r, dr_dtheta, found)

        type(geometry_t), intent(in) :: geom
        real(kind = rk), intent(in) :: cos_theta, sin_theta
        real(kind = rk), intent(out) :: r, dr_dtheta
        logical, intent(out) :: found

        integer(kind = ik), parameter :: ORDER(3) = [1_ik, 3_ik, 2_ik]
        real(kind = rk) :: a_sq, k, l, lead, roots(2), z_lo, z_hi, z, f_r, f_theta
        integer(kind = ik) :: j, seg, n_roots, m

        r = 0.0_rk
        dr_dtheta = 0.0_rk
        found = .false.

        do j = 1_ik, 3_ik
            seg = ORDER(j)
            a_sq = geom%a_sq(seg)
            k = geom%k(seg)
            l = geom%l(seg)
            select case (seg)
            case (1_ik)
                z_lo = -huge(1.0_rk)
                z_hi = geom%z_join(1) + RAY_JOIN_SLACK
            case (2_ik)
                z_lo = geom%z_join(2) - RAY_JOIN_SLACK
                z_hi = huge(1.0_rk)
            case default
                z_lo = geom%z_join(1) - RAY_JOIN_SLACK
                z_hi = geom%z_join(2) + RAY_JOIN_SLACK
            end select

            lead = sin_theta * sin_theta + k * cos_theta * cos_theta
            call solve_quadratic_s(lead, -2.0_rk * k * l * cos_theta, k * l * l - a_sq, &
                    n_roots, roots(1), roots(2))
            do m = 1_ik, min(n_roots, 2_ik)
                if (.not. (roots(m) > 0.0_rk)) cycle
                z = roots(m) * cos_theta
                if (z < z_lo .or. z > z_hi) cycle
                r = roots(m)
                f_r = 2.0_rk * r * lead - 2.0_rk * k * l * cos_theta
                f_theta = 2.0_rk * r * sin_theta * (r * cos_theta * (1.0_rk - k) + k * l)
                if (abs(f_r) > 0.0_rk) dr_dtheta = -f_theta / f_r
                found = .true.
                return
            end do
        end do

    end subroutine ray_radius_s

    !> Star-convexity margin from the origin: min over the shape of
    !! f(z) = rho - z drho/dz, closed form.
    !!
    !! On a segment rho rho' = -k (z - l), so
    !!   f = (rho^2 + k z (z - l))/rho = (a_sq + k l (z - l))/rho, and
    !!   f' = [k l rho^2 + (a_sq + k l (z - l)) k (z - l)]/rho^3 = k a_sq z/rho^3:
    !! the only interior critical point of any segment is z = 0, where
    !! f = rho(0). At the poles of a connected shape containing the origin f
    !! tends to +infinity. The minimum is therefore among f at the two joins
    !! (one segment each suffices, the joins are C1) and rho(0); f(0) is an
    !! attained value, so including it is safe where it is a maximum (a
    !! one-sheet middle). Exact, no sampling.
    !!
    !! @param[in] geom  Resolved, connected geometry (not separated)
    !! @return          The margin; the R(theta) outputs need >= STAR_CONVEXITY_MARGIN
    pure function star_margin_f(geom) result(margin)

        type(geometry_t), intent(in) :: geom
        real(kind = rk) :: margin

        real(kind = rk) :: rho, drho_dz

        call segment_rho_s(geom%a_sq(1), geom%k(1), geom%l(1), geom%z_join(1), rho, drho_dz)
        margin = rho - geom%z_join(1) * drho_dz
        call segment_rho_s(geom%a_sq(2), geom%k(2), geom%l(2), geom%z_join(2), rho, drho_dz)
        margin = min(margin, rho - geom%z_join(2) * drho_dz)
        call rho_at_z_s(geom, 0.0_rk, rho, drho_dz)
        margin = min(margin, rho)

    end function star_margin_f

    !===========================================================================
    ! MOLLER COORDINATES
    !===========================================================================

    !> Moller's grid coordinates of a resolved shape (R0 = 1 units).
    !!
    !! q2 = <r^2 P2(cos theta)> over the volume = 3/4 int (z^2 rho^2 - rho^4/4) dz.
    !! d: 2 a3 for a necked one-sheet middle, 0 when pinched or separated,
    !! otherwise the narrowest diameter of the middle body, which sits at a
    !! join (a spheroid arc is narrowest at an end, a flare is monotone).
    !! alpha_g = (a1^2 c1 - a2^2 c2)/(a1^2 c1 + a2^2 c2) (eq. 10).
    !! eps_i = 3 (1 - a_i/c_i)/(2 + a_i/c_i), inverting eqs. 6-7, a_i/c_i = sqrt(k_i).
    pure subroutine moller_coordinates_s(geom, q2, d, alpha_g, eps_f1, eps_f2)

        type(geometry_t), intent(in) :: geom
        real(kind = rk), intent(out) :: q2, d, alpha_g, eps_f1, eps_f2

        real(kind = rk) :: i0, i1, i2, i4, v1, v2, r1, r2, dummy

        call body_integrals_s(geom, i0, i1, i2, i4)
        q2 = 0.75_rk * (i2 - 0.25_rk * i4)

        if (geom%is_pinched .or. geom%is_separated) then
            d = 0.0_rk
        else if (geom%mid_kind == TQS_MID_ONE_SHEET .and. geom%waist_inside) then
            d = 2.0_rk * sqrt(geom%a_sq(3))
        else
            call rho_at_z_s(geom, geom%z_join(1), r1, dummy)
            call rho_at_z_s(geom, geom%z_join(2), r2, dummy)
            d = 2.0_rk * min(r1, r2)
        end if

        v1 = geom%a_sq(1) * geom%c(1)
        v2 = geom%a_sq(2) * geom%c(2)
        alpha_g = (v1 - v2) / (v1 + v2)

        eps_f1 = 3.0_rk * (1.0_rk - sqrt(geom%k(1))) / (2.0_rk + sqrt(geom%k(1)))
        eps_f2 = 3.0_rk * (1.0_rk - sqrt(geom%k(2))) / (2.0_rk + sqrt(geom%k(2)))

    end subroutine moller_coordinates_s

end module tqs_workers_mod
