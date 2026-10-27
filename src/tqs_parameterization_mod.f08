!> Public surface of the three-quadratic-surface (3QS) nuclear shape library.
!!
!! The library's front door: `cache_t`, the usage checks, the shape gates,
!! and every cached and one-shot routine. It does no floating-point
!! arithmetic on an output value; the code here checks, calls kernels,
!! applies gates and copies results. Every kernel lives in `tqs_resolve_mod`
!! (coordinates -> geometry) or `tqs_workers_mod` (arithmetic on a geometry).
!!
!! ## The shape
!!
!! Nix's three quadratic surfaces of revolution (Nucl. Phys. A 130 (1969)
!! 241, eq. 4): two end spheroids joined smoothly by a middle body that is a
!! spheroid, cylinder, hyperboloid of one sheet, cone or hyperboloid of two
!! sheets. `params = [sigma1, sigma2, sigma3, alpha2, alpha3]`, Nix's
!! coordinates, fixed length TQS_N_PARAMS = 5. The sphere is [0, 1, 1, 0, 0]
!! (and [0, s2, 1, 0, 0] for any s2). Everything is in R0 = 1 units with the
!! centre of mass at the origin; the library has no notion of A, Z or N.
!!
!! ## Two tiers
!!
!! 1. **One-shot** — `compute_*_standalone_s`: builds a local cache (the
!!    radius forms from the caller's thetas, the others from a one-node set),
!!    calls the cached routine and lets the cache go out of scope.
!! 2. **Read-only cache** — `cache_*_s` on a `cache_t` built once by
!!    `cache_init_s`. The cache holds only the primary theta set with its cos
!!    and sin; nothing derived from a shape parameter is ever stored.
!!
!! One-shot and cached calls return the same status and bitwise-identical
!! outputs: there is one pipeline, and the one-shot forms are wrappers.
!!
!! ## Thread and lifetime rules (normative)
!!
!! - A cache is immutable after `cache_init_s`; every compute takes it
!!   `intent(in)` and is `pure`. Any number of threads may compute on one
!!   cache concurrently.
!! - `cache_init_s` and `cache_free_s` on a cache must not race with any
!!   other call on it.
!! - No hidden state: the library holds no mutable module-level state.
!! - A `cache_t` has allocatable components only: assignment is a deep copy
!!   and scope exit releases it. `cache_free_s` is available, not mandatory.
!!
!! ## Checks, in order (normative)
!!
!! The first failing check sets the status. Usage codes depend on sizes,
!! handles and grids only, never on parameter values.
!!
!!   cache_init_s   max_params > 5 -> 1; < 5 -> 5; empty theta set or a
!!                  theta outside [0, pi] -> 3
!!   cached call    uninitialized cache -> 2; size(params) /= 5 -> 4; output
!!                  extent -> 3; at-thetas only: theta set -> 3
!!   one-shot call  size(params) > 5 -> 1; < 5 (empty included) -> 4; output
!!                  extent and theta set -> 3
!!
!! Then the value codes: 100 (coordinates), 101 (no geometry), then the
!! gates each output applies — 103 (separated), then 102 (pinched, or
!! star-convexity margin below STAR_CONVEXITY_MARGIN).
!!
!!   radius_grid, radius_and_derivative, at-thetas, shape   100 101 103 102
!!   rho_z_grid                                             100 101 103
!!   rho_z_grid_unchecked, neck, geometry, moller           100 101
!!   scission_sigma1                                        100
!!
!! On any nonzero status every output argument is zero-filled over its
!! entire extent, `found` and `defined` are .false., and a geometry record
!! is zero-filled with mid_kind = TQS_MID_NONE and every flag .false.
module tqs_parameterization_mod

    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use shape_core_mod, only: SHAPE_VALID, SHAPE_ERROR_TOO_MANY_PARAMS, &
            SHAPE_ERROR_CACHE_NOT_INITIALIZED, SHAPE_ERROR_INVALID_GRID, &
            SHAPE_ERROR_WRONG_PARAM_COUNT, SHAPE_ERROR_INVALID_INIT, &
            SHAPE_MAX_PARAMS
    use tqs_workers_mod, only: geometry_t, profile_s, neck_s, poles_s, &
            moller_coordinates_s, ray_radius_s, star_margin_f, theta_trig_s, &
            TQS_N_PARAMS, TQS_N_POINTS_FLOOR, &
            STAR_CONVEXITY_MARGIN, TQS_PINCH_TOL, TQS_DEGENERATE_TOL, &
            TQS_COINCIDENCE_TOL, TQS_MID_NONE, TQS_MID_SPHEROID, &
            TQS_MID_CYLINDER, TQS_MID_ONE_SHEET, TQS_MID_CONE, TQS_MID_TWO_SHEET, &
            TQS_ERROR_INVALID_COORDINATES, TQS_ERROR_NO_GEOMETRY, &
            TQS_ERROR_NOT_STAR_CONVEX, TQS_ERROR_SEPARATED
    use tqs_resolve_mod, only: resolve_s, alpha1_f, scission_sigma1_s

    implicit none

    private

    !---------------------------------------------------------------------------
    ! Read-only cache
    !---------------------------------------------------------------------------
    public :: cache_init_s
    public :: cache_free_s
    public :: cache_max_params_f
    public :: cache_n_thetas_f
    public :: cache_is_initialized_f

    !---------------------------------------------------------------------------
    ! Cached tier
    !---------------------------------------------------------------------------
    public :: cache_radius_grid_s
    public :: cache_radius_and_derivative_s
    public :: cache_radius_and_derivative_at_thetas_s
    public :: cache_shape_s
    public :: cache_rho_z_grid_s
    public :: cache_rho_z_grid_unchecked_s
    public :: cache_neck_s
    public :: cache_geometry_s
    public :: cache_moller_coordinates_s
    public :: cache_scission_sigma1_s

    !---------------------------------------------------------------------------
    ! One-shot tier
    !---------------------------------------------------------------------------
    public :: compute_radius_grid_standalone_s
    public :: compute_radius_and_derivative_standalone_s
    public :: compute_shape_standalone_s
    public :: compute_rho_z_grid_standalone_s
    public :: compute_rho_z_grid_unchecked_standalone_s
    public :: compute_neck_standalone_s
    public :: compute_geometry_standalone_s
    public :: compute_moller_coordinates_standalone_s
    public :: compute_scission_sigma1_standalone_s

    !---------------------------------------------------------------------------
    ! Constants, owned by the workers module and re-exported so the public
    ! surface is self-contained
    !---------------------------------------------------------------------------
    public :: TQS_N_PARAMS
    public :: TQS_N_POINTS_FLOOR
    public :: STAR_CONVEXITY_MARGIN
    public :: TQS_PINCH_TOL
    public :: TQS_DEGENERATE_TOL
    public :: TQS_COINCIDENCE_TOL
    public :: TQS_MID_NONE
    public :: TQS_MID_SPHEROID
    public :: TQS_MID_CYLINDER
    public :: TQS_MID_ONE_SHEET
    public :: TQS_MID_CONE
    public :: TQS_MID_TWO_SHEET
    public :: SHAPE_MAX_PARAMS

    !---------------------------------------------------------------------------
    ! Status codes
    !---------------------------------------------------------------------------
    !! Two disjoint ranges share one integer space. 0-99 are the shared
    !! shape_core codes, re-exported so a caller never has to use
    !! shape_core_mod directly; 100+ are tqs-specific. TQS_VALID is an alias
    !! of SHAPE_VALID, so success compares equal across libraries.
    public :: SHAPE_VALID
    public :: SHAPE_ERROR_TOO_MANY_PARAMS
    public :: SHAPE_ERROR_CACHE_NOT_INITIALIZED
    public :: SHAPE_ERROR_INVALID_GRID
    public :: SHAPE_ERROR_WRONG_PARAM_COUNT
    public :: SHAPE_ERROR_INVALID_INIT

    public :: TQS_ERROR_INVALID_COORDINATES
    public :: TQS_ERROR_NO_GEOMETRY
    public :: TQS_ERROR_NOT_STAR_CONVEX
    public :: TQS_ERROR_SEPARATED

    integer(kind = ik), parameter, public :: TQS_VALID = SHAPE_VALID

    !> Length of every string returned by status_message.
    integer(kind = ik), parameter, public :: STATUS_MESSAGE_LEN = 64_ik
    public :: status_message

    !> The one-node theta set of the theta-less one-shot forms.
    real(kind = rk), parameter :: ONE_THETA(1) = [0.0_rk]

    !> Read-only cache: the primary theta set with its cos and sin.
    !!
    !! Opaque: the components are private and may grow. Immutable after
    !! `cache_init_s`, shared freely across threads.
    type, public :: cache_t
        private
        logical :: initialized = .false.
        integer(kind = ik) :: max_params = 0_ik
        real(kind = rk), allocatable :: thetas(:)
        real(kind = rk), allocatable :: cos_t(:)
        real(kind = rk), allocatable :: sin_t(:)
    end type cache_t

    !> The resolved geometry as a diagnostic, in R0 = 1 units, centre of mass
    !! at the origin. Segment 1 is the south end spheroid, 2 the north one,
    !! 3 the middle body: rho^2 = a_sq(i) - k(i) (z - l(i))^2. c(1:2) are the
    !! end semi-axes along z, z_join(1:2) the joins, r_north and r_south the
    !! pole distances, alpha1 Nix's centre-of-mass coordinate. Succeeds for
    !! every valid geometry, separated and pinched included.
    type, public :: tqs_geometry_t
        real(kind = rk) :: u = 0.0_rk
        real(kind = rk) :: a_sq(3) = 0.0_rk
        real(kind = rk) :: k(3) = 0.0_rk
        real(kind = rk) :: c(2) = 0.0_rk
        real(kind = rk) :: l(3) = 0.0_rk
        real(kind = rk) :: z_join(2) = 0.0_rk
        real(kind = rk) :: r_north = 0.0_rk
        real(kind = rk) :: r_south = 0.0_rk
        real(kind = rk) :: alpha1 = 0.0_rk
        integer(kind = ik) :: mid_kind = TQS_MID_NONE
        logical :: waist_inside = .false.
        logical :: is_separated = .false.
        logical :: is_pinched = .false.
    end type tqs_geometry_t

contains

    !===========================================================================
    ! STATUS
    !===========================================================================

    !> Fixed diagnostic string for a status code, shared or tqs-specific.
    !!
    !! @param[in] code  Any integer; unrecognised values get a fallback message
    !! @return          Blank-padded message, never empty
    pure function status_message(code) result(msg)
        integer(kind = ik), intent(in) :: code
        character(len = STATUS_MESSAGE_LEN) :: msg
        select case (code)
        case (SHAPE_VALID);                        msg = 'valid'
        case (SHAPE_ERROR_TOO_MANY_PARAMS);        msg = 'too many parameters'
        case (SHAPE_ERROR_CACHE_NOT_INITIALIZED);  msg = 'cache not initialized'
        case (SHAPE_ERROR_INVALID_GRID);           msg = 'invalid grid: theta set, profile size or buffer extent'
        case (SHAPE_ERROR_WRONG_PARAM_COUNT);      msg = 'params length is not 5'
        case (SHAPE_ERROR_INVALID_INIT);           msg = 'invalid init arguments'
        case (TQS_ERROR_INVALID_COORDINATES);      msg = 'invalid coordinates: need k1, k2 > 0, |alpha2| < 2, sigma1 >= 0'
        case (TQS_ERROR_NO_GEOMETRY);              msg = 'no three-surface geometry for these coordinates'
        case (TQS_ERROR_NOT_STAR_CONVEX);          msg = 'shape not star-convex from the centre of mass'
        case (TQS_ERROR_SEPARATED);                msg = 'shape separated into two fragments'
        case default;                              msg = 'unknown status code'
        end select
    end function status_message

    !===========================================================================
    ! CACHE LIFECYCLE
    !===========================================================================

    !> Builds a read-only cache on a primary theta set.
    !!
    !! tqs is fixed-arity: `max_params` must be TQS_N_PARAMS; the argument
    !! stays for signature uniformity with the expansion libraries. Calling
    !! this on a live cache releases the old contents through `intent(out)`.
    !!
    !! @param[out] cache       Ready on success; uninitialized otherwise
    !! @param[in]  max_params  TQS_N_PARAMS (5)
    !! @param[in]  thetas      Primary polar angles in [0, pi], at least one, any order
    !! @param[out] status      SHAPE_VALID, SHAPE_ERROR_TOO_MANY_PARAMS (> 5),
    !!                         SHAPE_ERROR_INVALID_INIT (< 5) or SHAPE_ERROR_INVALID_GRID
    pure subroutine cache_init_s(cache, max_params, thetas, status)

        type(cache_t), intent(out) :: cache
        integer(kind = ik), intent(in) :: max_params
        real(kind = rk), intent(in) :: thetas(:)
        integer(kind = ik), intent(out) :: status

        integer(kind = ik) :: n
        integer :: alloc_stat

        if (max_params > TQS_N_PARAMS) then
            status = SHAPE_ERROR_TOO_MANY_PARAMS
            return
        end if

        if (max_params < TQS_N_PARAMS) then
            status = SHAPE_ERROR_INVALID_INIT
            return
        end if

        call check_thetas_s(thetas, status)
        if (status /= SHAPE_VALID) return

        n = size(thetas, kind = ik)
        allocate(cache%thetas(n), cache%cos_t(n), cache%sin_t(n), stat = alloc_stat)
        if (alloc_stat /= 0) then
            status = SHAPE_ERROR_INVALID_GRID
            return
        end if

        cache%thetas = thetas
        call theta_trig_s(thetas, cache%cos_t, cache%sin_t)
        cache%max_params = max_params
        cache%initialized = .true.

    end subroutine cache_init_s

    !> Releases the cache. Infallible, and safe on an already-free instance.
    pure subroutine cache_free_s(cache)

        type(cache_t), intent(inout) :: cache

        if (allocated(cache%thetas)) deallocate(cache%thetas)
        if (allocated(cache%cos_t)) deallocate(cache%cos_t)
        if (allocated(cache%sin_t)) deallocate(cache%sin_t)
        cache%max_params = 0_ik
        cache%initialized = .false.

    end subroutine cache_free_s

    !> The cache's max_params; 0 when uninitialized.
    pure function cache_max_params_f(cache) result(max_params)
        type(cache_t), intent(in) :: cache
        integer(kind = ik) :: max_params
        max_params = cache%max_params
    end function cache_max_params_f

    !> Number of primary thetas of the cache; 0 when uninitialized.
    pure function cache_n_thetas_f(cache) result(n_thetas)
        type(cache_t), intent(in) :: cache
        integer(kind = ik) :: n_thetas
        n_thetas = 0_ik
        if (cache%initialized) n_thetas = size(cache%thetas, kind = ik)
    end function cache_n_thetas_f

    !> .true. after a successful `cache_init_s`, .false. after `cache_free_s`.
    pure function cache_is_initialized_f(cache) result(is_initialized)
        type(cache_t), intent(in) :: cache
        logical :: is_initialized
        is_initialized = cache%initialized
    end function cache_is_initialized_f

    !===========================================================================
    ! CACHED TIER
    !===========================================================================

    !> R(theta) at the cache's primary thetas.
    !!
    !! @param[in]  cache   Initialized cache
    !! @param[in]  params  [sigma1, sigma2, sigma3, alpha2, alpha3]
    !! @param[out] radii   R(theta), cache_n_thetas_f(cache) long
    !! @param[out] status  0, 2, 4, 3, 100, 101, 103 or 102
    pure subroutine cache_radius_grid_s(cache, params, radii, status)

        type(cache_t), intent(in) :: cache
        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(out) :: radii(:)
        integer(kind = ik), intent(out) :: status

        real(kind = rk), allocatable :: dr_sink(:)
        integer :: alloc_stat

        radii = 0.0_rk

        call check_call_s(cache, params, status)
        if (status /= SHAPE_VALID) return

        if (size(radii, kind = ik) /= size(cache%thetas, kind = ik)) then
            status = SHAPE_ERROR_INVALID_GRID
            return
        end if

        ! The ray solve always produces the derivative; this output does not
        ! report it, so it lands in a write-only sink.
        allocate(dr_sink(size(radii)), stat = alloc_stat)
        if (alloc_stat /= 0) then
            status = SHAPE_ERROR_INVALID_GRID
            return
        end if

        call resolved_radii_s(params, cache%cos_t, cache%sin_t, radii, dr_sink, status)

    end subroutine cache_radius_grid_s

    !> R(theta) and dR/dtheta at the cache's primary thetas.
    !!
    !! @param[in]  cache      Initialized cache
    !! @param[in]  params     [sigma1, sigma2, sigma3, alpha2, alpha3]
    !! @param[out] radii      R(theta), cache_n_thetas_f(cache) long
    !! @param[out] dr_dtheta  dR/dtheta, same length
    !! @param[out] status     0, 2, 4, 3, 100, 101, 103 or 102
    pure subroutine cache_radius_and_derivative_s(cache, params, radii, dr_dtheta, status)

        type(cache_t), intent(in) :: cache
        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(out) :: radii(:)
        real(kind = rk), intent(out) :: dr_dtheta(:)
        integer(kind = ik), intent(out) :: status

        call zero_fill_pair_s(radii, dr_dtheta)

        call check_call_s(cache, params, status)
        if (status /= SHAPE_VALID) return

        if (size(radii, kind = ik) /= size(cache%thetas, kind = ik) &
                .or. size(dr_dtheta, kind = ik) /= size(cache%thetas, kind = ik)) then
            status = SHAPE_ERROR_INVALID_GRID
            return
        end if

        call resolved_radii_s(params, cache%cos_t, cache%sin_t, radii, dr_dtheta, status)

    end subroutine cache_radius_and_derivative_s

    !> R(theta) and dR/dtheta at caller-supplied thetas (the WMMM hot path).
    !!
    !! Given the cache's own thetas it returns the same bits as
    !! `cache_radius_and_derivative_s`: the same trig kernel, the same solve.
    !!
    !! @param[in]  cache      Initialized cache
    !! @param[in]  params     [sigma1, sigma2, sigma3, alpha2, alpha3]
    !! @param[in]  thetas     Polar angles in [0, pi], at least one, any order
    !! @param[out] radii      R(theta), size(thetas) long
    !! @param[out] dr_dtheta  dR/dtheta, size(thetas) long
    !! @param[out] status     0, 2, 4, 3, 100, 101, 103 or 102
    pure subroutine cache_radius_and_derivative_at_thetas_s(cache, params, thetas, &
            radii, dr_dtheta, status)

        type(cache_t), intent(in) :: cache
        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(in) :: thetas(:)
        real(kind = rk), intent(out) :: radii(:)
        real(kind = rk), intent(out) :: dr_dtheta(:)
        integer(kind = ik), intent(out) :: status

        real(kind = rk), allocatable :: cos_t(:), sin_t(:)
        integer(kind = ik) :: n
        integer :: alloc_stat

        call zero_fill_pair_s(radii, dr_dtheta)

        call check_call_s(cache, params, status)
        if (status /= SHAPE_VALID) return

        n = size(thetas, kind = ik)
        if (size(radii, kind = ik) /= n .or. size(dr_dtheta, kind = ik) /= n) then
            status = SHAPE_ERROR_INVALID_GRID
            return
        end if

        call check_thetas_s(thetas, status)
        if (status /= SHAPE_VALID) return

        allocate(cos_t(n), sin_t(n), stat = alloc_stat)
        if (alloc_stat /= 0) then
            status = SHAPE_ERROR_INVALID_GRID
            return
        end if
        call theta_trig_s(thetas, cos_t, sin_t)

        call resolved_radii_s(params, cos_t, sin_t, radii, dr_dtheta, status)

    end subroutine cache_radius_and_derivative_at_thetas_s

    !> Pole radii in the centre-of-mass frame, with the R(theta) verdict.
    !!
    !! @param[in]  cache    Initialized cache
    !! @param[in]  params   [sigma1, sigma2, sigma3, alpha2, alpha3]
    !! @param[out] r_north  R(0)
    !! @param[out] r_south  R(pi)
    !! @param[out] status   0, 2, 4, 100, 101, 103 or 102
    pure subroutine cache_shape_s(cache, params, r_north, r_south, status)

        type(cache_t), intent(in) :: cache
        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(out) :: r_north, r_south
        integer(kind = ik), intent(out) :: status

        type(geometry_t) :: geom

        r_north = 0.0_rk
        r_south = 0.0_rk

        call check_call_s(cache, params, status)
        if (status /= SHAPE_VALID) return

        call resolve_s(params, geom, status)
        if (status /= SHAPE_VALID) return

        call gate_rtheta_s(geom, status)
        if (status /= SHAPE_VALID) return

        call poles_s(geom, r_north, r_south)

    end subroutine cache_shape_s

    !> rho(z) and drho/dz on a uniform grid over [z_south, z_north], both poles
    !! included exactly; (0, 0) at the poles and at a pinch point. Checked:
    !! a separated shape is 103.
    !!
    !! @param[in]  cache    Initialized cache
    !! @param[in]  params   [sigma1, sigma2, sigma3, alpha2, alpha3]
    !! @param[out] z        Grid nodes; size >= TQS_N_POINTS_FLOOR
    !! @param[out] rho      rho(z), size(z) long
    !! @param[out] drho_dz  drho/dz, size(z) long
    !! @param[out] status   0, 2, 4, 3, 100, 101 or 103
    pure subroutine cache_rho_z_grid_s(cache, params, z, rho, drho_dz, status)

        type(cache_t), intent(in) :: cache
        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(out) :: z(:), rho(:), drho_dz(:)
        integer(kind = ik), intent(out) :: status

        call rho_z_grid_core_s(cache, params, .true., z, rho, drho_dz, status)

    end subroutine cache_rho_z_grid_s

    !> As `cache_rho_z_grid_s`, never 103: a separated shape is drawn with
    !! (0, 0) at every node inside the void.
    pure subroutine cache_rho_z_grid_unchecked_s(cache, params, z, rho, drho_dz, status)

        type(cache_t), intent(in) :: cache
        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(out) :: z(:), rho(:), drho_dz(:)
        integer(kind = ik), intent(out) :: status

        call rho_z_grid_core_s(cache, params, .false., z, rho, drho_dz, status)

    end subroutine cache_rho_z_grid_unchecked_s

    !> The neck of a middle body with sigma2 < 0 whose waist lies strictly
    !! inside it: z_neck = l3, rho_neck = a3, or 0 for a pinched or separated
    !! shape (the split plane). found = .false. otherwise.
    !!
    !! @param[out] status  0, 2, 4, 100 or 101
    pure subroutine cache_neck_s(cache, params, z_neck, rho_neck, found, status)

        type(cache_t), intent(in) :: cache
        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(out) :: z_neck, rho_neck
        logical, intent(out) :: found
        integer(kind = ik), intent(out) :: status

        type(geometry_t) :: geom

        z_neck = 0.0_rk
        rho_neck = 0.0_rk
        found = .false.

        call check_call_s(cache, params, status)
        if (status /= SHAPE_VALID) return

        call resolve_s(params, geom, status)
        if (status /= SHAPE_VALID) return

        call neck_s(geom, z_neck, rho_neck, found)

    end subroutine cache_neck_s

    !> The resolved geometry; succeeds for every valid geometry, separated and
    !! pinched included, so it is the source of poles and flags when the
    !! R(theta) outputs report 102 or 103.
    !!
    !! @param[out] geometry  The record; zero-filled with mid_kind = TQS_MID_NONE on failure
    !! @param[out] status    0, 2, 4, 100 or 101
    pure subroutine cache_geometry_s(cache, params, geometry, status)

        type(cache_t), intent(in) :: cache
        real(kind = rk), intent(in) :: params(:)
        type(tqs_geometry_t), intent(out) :: geometry
        integer(kind = ik), intent(out) :: status

        type(geometry_t) :: geom

        call check_call_s(cache, params, status)
        if (status /= SHAPE_VALID) return

        call resolve_s(params, geom, status)
        if (status /= SHAPE_VALID) return

        geometry%u = geom%u
        geometry%a_sq = geom%a_sq
        geometry%k = geom%k
        geometry%c = geom%c
        geometry%l = geom%l
        geometry%z_join = geom%z_join
        call poles_s(geom, geometry%r_north, geometry%r_south)
        geometry%alpha1 = alpha1_f(geom)
        geometry%mid_kind = geom%mid_kind
        geometry%waist_inside = geom%waist_inside
        geometry%is_separated = geom%is_separated
        geometry%is_pinched = geom%is_pinched

    end subroutine cache_geometry_s

    !> Moller's grid coordinates of the resolved shape, as diagnostics: q2 =
    !! <r^2 P2> (R0 = 1), d, alpha_g, eps_f1, eps_f2. d = 0 for pinched and
    !! separated shapes.
    !!
    !! @param[out] status  0, 2, 4, 100 or 101
    pure subroutine cache_moller_coordinates_s(cache, params, q2, d, alpha_g, eps_f1, &
            eps_f2, status)

        type(cache_t), intent(in) :: cache
        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(out) :: q2, d, alpha_g, eps_f1, eps_f2
        integer(kind = ik), intent(out) :: status

        type(geometry_t) :: geom

        q2 = 0.0_rk
        d = 0.0_rk
        alpha_g = 0.0_rk
        eps_f1 = 0.0_rk
        eps_f2 = 0.0_rk

        call check_call_s(cache, params, status)
        if (status /= SHAPE_VALID) return

        call resolve_s(params, geom, status)
        if (status /= SHAPE_VALID) return

        call moller_coordinates_s(geom, q2, d, alpha_g, eps_f1, eps_f2)

    end subroutine cache_moller_coordinates_s

    !> The scission map of (sigma2, sigma3, alpha2, alpha3): Moller's
    !! touching-spheroids sigma1_t always, and the a3^2 = 0 locus sigma1_star
    !! with `defined` = .false. (and 0) for sigma2 >= 0. sigma1 is validated
    !! but does not enter.
    !!
    !! @param[out] status  0, 2, 4 or 100
    pure subroutine cache_scission_sigma1_s(cache, params, sigma1_t, sigma1_star, defined, &
            status)

        type(cache_t), intent(in) :: cache
        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(out) :: sigma1_t, sigma1_star
        logical, intent(out) :: defined
        integer(kind = ik), intent(out) :: status

        sigma1_t = 0.0_rk
        sigma1_star = 0.0_rk
        defined = .false.

        call check_call_s(cache, params, status)
        if (status /= SHAPE_VALID) return

        call scission_sigma1_s(params, sigma1_t, sigma1_star, defined, status)

    end subroutine cache_scission_sigma1_s

    !===========================================================================
    ! ONE-SHOT TIER
    !===========================================================================

    !> R(theta) at caller-supplied thetas, one-shot.
    !!
    !! @param[in]  params  [sigma1, sigma2, sigma3, alpha2, alpha3]
    !! @param[in]  thetas  Polar angles in [0, pi], at least one
    !! @param[out] radii   R(theta), size(thetas) long
    !! @param[out] status  As the cached form, with 1 for > 5 params
    pure subroutine compute_radius_grid_standalone_s(params, thetas, radii, status)

        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(in) :: thetas(:)
        real(kind = rk), intent(out) :: radii(:)
        integer(kind = ik), intent(out) :: status

        type(cache_t) :: cache

        radii = 0.0_rk

        call check_one_shot_s(params, status)
        if (status /= SHAPE_VALID) return

        if (size(radii, kind = ik) /= size(thetas, kind = ik)) then
            status = SHAPE_ERROR_INVALID_GRID
            return
        end if

        call cache_init_s(cache, TQS_N_PARAMS, thetas, status)
        if (status /= SHAPE_VALID) return

        call cache_radius_grid_s(cache, params, radii, status)

    end subroutine compute_radius_grid_standalone_s

    !> R(theta) and dR/dtheta at caller-supplied thetas, one-shot.
    pure subroutine compute_radius_and_derivative_standalone_s(params, thetas, radii, &
            dr_dtheta, status)

        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(in) :: thetas(:)
        real(kind = rk), intent(out) :: radii(:)
        real(kind = rk), intent(out) :: dr_dtheta(:)
        integer(kind = ik), intent(out) :: status

        type(cache_t) :: cache

        call zero_fill_pair_s(radii, dr_dtheta)

        call check_one_shot_s(params, status)
        if (status /= SHAPE_VALID) return

        if (size(radii, kind = ik) /= size(thetas, kind = ik) &
                .or. size(dr_dtheta, kind = ik) /= size(thetas, kind = ik)) then
            status = SHAPE_ERROR_INVALID_GRID
            return
        end if

        call cache_init_s(cache, TQS_N_PARAMS, thetas, status)
        if (status /= SHAPE_VALID) return

        call cache_radius_and_derivative_s(cache, params, radii, dr_dtheta, status)

    end subroutine compute_radius_and_derivative_standalone_s

    !> Pole radii with the R(theta) verdict, one-shot.
    pure subroutine compute_shape_standalone_s(params, r_north, r_south, status)

        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(out) :: r_north, r_south
        integer(kind = ik), intent(out) :: status

        type(cache_t) :: cache

        r_north = 0.0_rk
        r_south = 0.0_rk

        call one_shot_cache_s(params, cache, status)
        if (status /= SHAPE_VALID) return

        call cache_shape_s(cache, params, r_north, r_south, status)

    end subroutine compute_shape_standalone_s

    !> Checked rho(z) profile, one-shot.
    pure subroutine compute_rho_z_grid_standalone_s(params, z, rho, drho_dz, status)

        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(out) :: z(:), rho(:), drho_dz(:)
        integer(kind = ik), intent(out) :: status

        type(cache_t) :: cache

        call zero_fill_grid_s(z, rho, drho_dz)

        call one_shot_cache_s(params, cache, status)
        if (status /= SHAPE_VALID) return

        call cache_rho_z_grid_s(cache, params, z, rho, drho_dz, status)

    end subroutine compute_rho_z_grid_standalone_s

    !> Unchecked rho(z) profile, one-shot.
    pure subroutine compute_rho_z_grid_unchecked_standalone_s(params, z, rho, drho_dz, status)

        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(out) :: z(:), rho(:), drho_dz(:)
        integer(kind = ik), intent(out) :: status

        type(cache_t) :: cache

        call zero_fill_grid_s(z, rho, drho_dz)

        call one_shot_cache_s(params, cache, status)
        if (status /= SHAPE_VALID) return

        call cache_rho_z_grid_unchecked_s(cache, params, z, rho, drho_dz, status)

    end subroutine compute_rho_z_grid_unchecked_standalone_s

    !> Neck, one-shot.
    pure subroutine compute_neck_standalone_s(params, z_neck, rho_neck, found, status)

        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(out) :: z_neck, rho_neck
        logical, intent(out) :: found
        integer(kind = ik), intent(out) :: status

        type(cache_t) :: cache

        z_neck = 0.0_rk
        rho_neck = 0.0_rk
        found = .false.

        call one_shot_cache_s(params, cache, status)
        if (status /= SHAPE_VALID) return

        call cache_neck_s(cache, params, z_neck, rho_neck, found, status)

    end subroutine compute_neck_standalone_s

    !> Resolved geometry, one-shot.
    pure subroutine compute_geometry_standalone_s(params, geometry, status)

        real(kind = rk), intent(in) :: params(:)
        type(tqs_geometry_t), intent(out) :: geometry
        integer(kind = ik), intent(out) :: status

        type(cache_t) :: cache

        call one_shot_cache_s(params, cache, status)
        if (status /= SHAPE_VALID) return

        call cache_geometry_s(cache, params, geometry, status)

    end subroutine compute_geometry_standalone_s

    !> Moller coordinates, one-shot.
    pure subroutine compute_moller_coordinates_standalone_s(params, q2, d, alpha_g, eps_f1, &
            eps_f2, status)

        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(out) :: q2, d, alpha_g, eps_f1, eps_f2
        integer(kind = ik), intent(out) :: status

        type(cache_t) :: cache

        q2 = 0.0_rk
        d = 0.0_rk
        alpha_g = 0.0_rk
        eps_f1 = 0.0_rk
        eps_f2 = 0.0_rk

        call one_shot_cache_s(params, cache, status)
        if (status /= SHAPE_VALID) return

        call cache_moller_coordinates_s(cache, params, q2, d, alpha_g, eps_f1, eps_f2, status)

    end subroutine compute_moller_coordinates_standalone_s

    !> Scission map, one-shot.
    pure subroutine compute_scission_sigma1_standalone_s(params, sigma1_t, sigma1_star, &
            defined, status)

        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(out) :: sigma1_t, sigma1_star
        logical, intent(out) :: defined
        integer(kind = ik), intent(out) :: status

        type(cache_t) :: cache

        sigma1_t = 0.0_rk
        sigma1_star = 0.0_rk
        defined = .false.

        call one_shot_cache_s(params, cache, status)
        if (status /= SHAPE_VALID) return

        call cache_scission_sigma1_s(cache, params, sigma1_t, sigma1_star, defined, status)

    end subroutine compute_scission_sigma1_standalone_s

    !===========================================================================
    ! CHECKS, PIPELINE, GATES
    !===========================================================================

    !> The two usage checks every cached compute starts with: an uninitialized
    !! cache (2), then a parameter count other than max_params (4).
    pure subroutine check_call_s(cache, params, status)

        type(cache_t), intent(in) :: cache
        real(kind = rk), intent(in) :: params(:)
        integer(kind = ik), intent(out) :: status

        if (.not. cache%initialized) then
            status = SHAPE_ERROR_CACHE_NOT_INITIALIZED
            return
        end if

        if (size(params, kind = ik) /= cache%max_params) then
            status = SHAPE_ERROR_WRONG_PARAM_COUNT
            return
        end if

        status = SHAPE_VALID

    end subroutine check_call_s

    !> The parameter-count check of every one-shot form: longer than
    !! TQS_N_PARAMS (1), shorter, empty included (4). Never padded, never trimmed.
    pure subroutine check_one_shot_s(params, status)

        real(kind = rk), intent(in) :: params(:)
        integer(kind = ik), intent(out) :: status

        if (size(params, kind = ik) > TQS_N_PARAMS) then
            status = SHAPE_ERROR_TOO_MANY_PARAMS
        else if (size(params, kind = ik) < TQS_N_PARAMS) then
            status = SHAPE_ERROR_WRONG_PARAM_COUNT
        else
            status = SHAPE_VALID
        end if

    end subroutine check_one_shot_s

    !> The count check and the one-node cache of the theta-less one-shot forms.
    pure subroutine one_shot_cache_s(params, cache, status)

        real(kind = rk), intent(in) :: params(:)
        type(cache_t), intent(out) :: cache
        integer(kind = ik), intent(out) :: status

        call check_one_shot_s(params, status)
        if (status /= SHAPE_VALID) return

        call cache_init_s(cache, TQS_N_PARAMS, ONE_THETA, status)

    end subroutine one_shot_cache_s

    !> A theta set: at least one node, every node in [0, pi].
    pure subroutine check_thetas_s(thetas, status)

        real(kind = rk), intent(in) :: thetas(:)
        integer(kind = ik), intent(out) :: status

        integer(kind = ik) :: i

        status = SHAPE_VALID

        if (size(thetas, kind = ik) < 1_ik) then
            status = SHAPE_ERROR_INVALID_GRID
            return
        end if

        do i = 1_ik, size(thetas, kind = ik)
            if (thetas(i) < 0.0_rk .or. thetas(i) > PI_C) then
                status = SHAPE_ERROR_INVALID_GRID
                return
            end if
        end do

    end subroutine check_thetas_s

    !> The R(theta) verdict on a resolved geometry: separated (103), then
    !! pinched (102), then a star-convexity margin below STAR_CONVEXITY_MARGIN (102).
    pure subroutine gate_rtheta_s(geom, status)

        type(geometry_t), intent(in) :: geom
        integer(kind = ik), intent(out) :: status

        status = SHAPE_VALID
        if (geom%is_separated) then
            status = TQS_ERROR_SEPARATED
        else if (geom%is_pinched) then
            status = TQS_ERROR_NOT_STAR_CONVEX
        else if (star_margin_f(geom) < STAR_CONVEXITY_MARGIN) then
            status = TQS_ERROR_NOT_STAR_CONVEX
        end if

    end subroutine gate_rtheta_s

    !> The shared body of every R(theta) output: resolve, gate, one ray per theta.
    !! On a nonzero status both outputs are zero-filled.
    pure subroutine resolved_radii_s(params, cos_t, sin_t, radii, dr_dtheta, status)

        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(in) :: cos_t(:), sin_t(:)
        real(kind = rk), intent(inout) :: radii(:), dr_dtheta(:)
        integer(kind = ik), intent(out) :: status

        type(geometry_t) :: geom
        integer(kind = ik) :: i
        logical :: found

        call resolve_s(params, geom, status)
        if (status /= SHAPE_VALID) return

        call gate_rtheta_s(geom, status)
        if (status /= SHAPE_VALID) return

        do i = 1_ik, size(cos_t, kind = ik)
            call ray_radius_s(geom, cos_t(i), sin_t(i), radii(i), dr_dtheta(i), found)
            ! Unreachable past the star gate; kept so a miss is a status, not a zero radius.
            if (.not. found) then
                status = TQS_ERROR_NOT_STAR_CONVEX
                call zero_fill_pair_s(radii, dr_dtheta)
                return
            end if
        end do

    end subroutine resolved_radii_s

    !> The shared body of the two cached profile outputs.
    pure subroutine rho_z_grid_core_s(cache, params, checked, z, rho, drho_dz, status)

        type(cache_t), intent(in) :: cache
        real(kind = rk), intent(in) :: params(:)
        logical, intent(in) :: checked
        real(kind = rk), intent(out) :: z(:), rho(:), drho_dz(:)
        integer(kind = ik), intent(out) :: status

        type(geometry_t) :: geom
        integer(kind = ik) :: n

        call zero_fill_grid_s(z, rho, drho_dz)

        call check_call_s(cache, params, status)
        if (status /= SHAPE_VALID) return

        n = size(z, kind = ik)
        if (n < TQS_N_POINTS_FLOOR .or. size(rho, kind = ik) /= n &
                .or. size(drho_dz, kind = ik) /= n) then
            status = SHAPE_ERROR_INVALID_GRID
            return
        end if

        call resolve_s(params, geom, status)
        if (status /= SHAPE_VALID) return

        if (checked .and. geom%is_separated) then
            status = TQS_ERROR_SEPARATED
            return
        end if

        call profile_s(geom, z, rho, drho_dz)

    end subroutine rho_z_grid_core_s

    !> Zeroes the full actual extent of a radius/derivative output pair.
    pure subroutine zero_fill_pair_s(a, b)

        real(kind = rk), intent(out) :: a(:)
        real(kind = rk), intent(out) :: b(:)

        a = 0.0_rk
        b = 0.0_rk

    end subroutine zero_fill_pair_s

    !> Zeroes the full actual extent of the profile outputs.
    pure subroutine zero_fill_grid_s(z, rho, drho_dz)

        real(kind = rk), intent(out) :: z(:)
        real(kind = rk), intent(out) :: rho(:)
        real(kind = rk), intent(out) :: drho_dz(:)

        z = 0.0_rk
        rho = 0.0_rk
        drho_dz = 0.0_rk

    end subroutine zero_fill_grid_s

end module tqs_parameterization_mod
