!> Status-code tests: tqs codes are >= 100, shared shape_core codes re-export,
!! every code maps to a non-empty message, and the fixed constants hold.
program tqs_param_status_test
    use precision_utilities_mod, only: ik, rk
    use tqs_parameterization_mod, only: TQS_VALID, TQS_ERROR_INVALID_COORDINATES, &
            TQS_ERROR_NO_GEOMETRY, TQS_ERROR_NOT_STAR_CONVEX, TQS_ERROR_SEPARATED, &
            TQS_N_PARAMS, TQS_N_POINTS_FLOOR, TQS_MID_SPHEROID, TQS_MID_CYLINDER, &
            TQS_MID_ONE_SHEET, TQS_MID_CONE, TQS_MID_TWO_SHEET, STAR_CONVEXITY_MARGIN, &
            SHAPE_VALID, SHAPE_ERROR_TOO_MANY_PARAMS, status_message
    use test_utils_mod, only: assert_int_eq, assert_true, assert_abs_close, test_summary
    implicit none

    call assert_int_eq(int(TQS_VALID), int(SHAPE_VALID), 'TQS_VALID aliases SHAPE_VALID')
    call assert_int_eq(int(TQS_ERROR_INVALID_COORDINATES), 100, 'invalid_coordinates = 100')
    call assert_int_eq(int(TQS_ERROR_NO_GEOMETRY), 101, 'no_geometry = 101')
    call assert_int_eq(int(TQS_ERROR_NOT_STAR_CONVEX), 102, 'not_star_convex = 102')
    call assert_int_eq(int(TQS_ERROR_SEPARATED), 103, 'separated = 103')
    call assert_int_eq(int(TQS_N_PARAMS), 5, 'five coordinates')
    call assert_int_eq(int(TQS_N_POINTS_FLOOR), 3, 'profile floor')
    call assert_int_eq(int(TQS_MID_SPHEROID), 1, 'mid spheroid = 1')
    call assert_int_eq(int(TQS_MID_CYLINDER), 2, 'mid cylinder = 2')
    call assert_int_eq(int(TQS_MID_ONE_SHEET), 3, 'mid one_sheet = 3')
    call assert_int_eq(int(TQS_MID_CONE), 4, 'mid cone = 4')
    call assert_int_eq(int(TQS_MID_TWO_SHEET), 5, 'mid two_sheet = 5')
    call assert_abs_close(STAR_CONVEXITY_MARGIN, 1.0e-2_rk, 0.0_rk, 'margin as fos')
    call assert_true(len_trim(status_message(101_ik)) > 0, 'message for 101 nonempty')
    call assert_true(len_trim(status_message(SHAPE_ERROR_TOO_MANY_PARAMS)) > 0, 'shared code has message')
    call assert_true(len_trim(status_message(999_ik)) > 0, 'unknown code has fallback message')
    call test_summary()
end program tqs_param_status_test
