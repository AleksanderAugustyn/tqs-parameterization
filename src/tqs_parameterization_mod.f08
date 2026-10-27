!> Public surface of the three-quadratic-surface (3QS) nuclear shape library.
!!
!! The library's front door: status codes, the re-exported constants, and
!! (from the front-module task on) the two tiers. It does no floating-point
!! arithmetic on an output value; every kernel lives in `tqs_workers_mod` or
!! `tqs_resolve_mod`.
module tqs_parameterization_mod

    use precision_utilities_mod, only: ik
    use shape_core_mod, only: SHAPE_VALID, SHAPE_ERROR_TOO_MANY_PARAMS, &
            SHAPE_ERROR_CACHE_NOT_INITIALIZED, SHAPE_ERROR_INVALID_GRID, &
            SHAPE_ERROR_WRONG_PARAM_COUNT, SHAPE_ERROR_INVALID_INIT, &
            SHAPE_MAX_PARAMS
    use tqs_workers_mod, only: TQS_N_PARAMS, TQS_N_POINTS_FLOOR, &
            STAR_CONVEXITY_MARGIN, TQS_PINCH_TOL, TQS_DEGENERATE_TOL, &
            TQS_COINCIDENCE_TOL, TQS_MID_NONE, TQS_MID_SPHEROID, TQS_MID_CYLINDER, &
            TQS_MID_ONE_SHEET, TQS_MID_CONE, TQS_MID_TWO_SHEET, &
            TQS_ERROR_INVALID_COORDINATES, TQS_ERROR_NO_GEOMETRY, &
            TQS_ERROR_NOT_STAR_CONVEX, TQS_ERROR_SEPARATED

    implicit none

    private

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

end module tqs_parameterization_mod
