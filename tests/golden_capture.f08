!> Prints golden values for the regression shapes as paste-ready Python.
!!
!! ## Recapture procedure
!!
!! The goldens in `python/tests/test_api.py` are comparisons against this
!! library at GOLDEN_ATOL = 1e-13. They must be regenerated whenever the
!! arithmetic behind the C API changes deliberately — a new resolve, a new
!! integral form, a new ray solve — and NOT otherwise: a golden that drifts
!! without a deliberate recapture is a regression, not a refresh.
!!
!!   1. Build the capture program (Debug):
!!        `cmake --build build/debug --target tqs_param_golden_capture`
!!   2. Run it and keep the output:
!!        `./build/debug/tqs_param_golden_capture`
!!      Every line but the closing `Tests:` summary is valid Python. A
!!      non-zero exit or a FAIL line means a golden shape stopped being
!!      accepted — investigate, do not paste.
!!   3. Replace the generated block of `python/tests/test_api.py` (between the
!!      `--- golden block` markers) with the printed lines.
!!   4. Run the Python suite against THIS tree, not an installed package:
!!        `TQS_PARAM_LIB=$PWD/build/release/libtqs_parameterization.so \
!!           PYTHONPATH=$PWD/python python -m pytest python/tests -q -p no:cacheprovider`
!!      Both variables are required: a bare `pytest` imports whatever
!!      `tqs_parameterization` is installed and loads whatever library that
!!      resolves. The suite must pass with GOLDEN_ATOL UNCHANGED; loosening it
!!      to make a recapture fit is never the fix.
!!
!! ## The tolerance budget
!!
!! Capture is Debug, but the goldens are compared against local Release
!! (-O3 -ffast-math) and the manylinux wheel (GCC 10.2, nehalem) as well. The
!! library is closed form — no iteration that self-corrects to a fixed point —
!! so every value inherits the FMA/contraction differences between those
!! builds directly: a few ulp of values of order unity. GOLDEN_ATOL = 1e-13
!! leaves two orders of magnitude for that gap; an algorithm change that eats
!! it must be measured against the shipped configuration, never answered by a
!! looser tolerance.
!!
!! Every value is captured through the SAME C entries the Python bindings
!! call, so the goldens include the C-level marshalling.
program golden_capture

    use precision_utilities_mod, only: ik, rk
    use c_bindings_mod, only: ik_c, rk_c
    use tqs_parameterization_c_api_mod, only: tqs_param_geometry_t, tqs_param_geometry, &
            tqs_param_neck, tqs_param_moller_coordinates, tqs_param_scission_sigma1, &
            tqs_param_radius_and_derivative, tqs_param_rho_z_grid
    use tqs_parameterization_mod, only: TQS_VALID
    use test_utils_mod, only: assert_int_eq, test_summary

    implicit none

    integer(kind = ik), parameter :: N_Z = 181_ik
    integer(kind = ik), parameter :: IDX(5) = [1_ik, 46_ik, 91_ik, 136_ik, 181_ik]

    call capture('SYM_NECKED', [1.6_rk, -0.6_rk, 1.0_rk, 0.0_rk, 0.0_rk])
    call capture('ASYM_NECKED', [1.6_rk, -0.6_rk, 1.0_rk, 0.3_rk, 0.2_rk])
    call capture('PEAR', [1.0_rk, -1.0_rk, 1.0_rk, -1.4_rk, 0.0_rk])
    call capture('SPHEROID_MID', [0.8_rk, 0.5_rk, 1.0_rk, 0.2_rk, 0.1_rk])
    call test_summary()

contains

    subroutine print_list(name, suffix, values)
        character(len = *), intent(in) :: name, suffix
        real(kind = rk_c), intent(in) :: values(:)
        integer(kind = ik) :: i
        write(*, '(A,A,A)') name, suffix, ' = ['
        do i = 1_ik, size(values, kind = ik)
            write(*, '(A,ES24.16E3,A)') '    ', values(i), ','
        end do
        write(*, '(A)') ']'
    end subroutine print_list

    subroutine capture(name, params)
        character(len = *), intent(in) :: name
        real(kind = rk), intent(in) :: params(5)

        real(kind = rk_c) :: p(5), geo(8), nk(2), mo(5), sc(2), z(N_Z), rho(N_Z), drho(N_Z)
        real(kind = rk_c) :: eval_thetas(3), eval_r(3), eval_dr(3), picked(size(IDX))
        type(tqs_param_geometry_t) :: g
        integer(kind = ik_c) :: code, found, defined

        p = real(params, rk_c)

        call tqs_param_geometry(p, 5_ik_c, g, code)
        call assert_int_eq(int(code, ik), TQS_VALID, name // ': geometry valid')
        if (int(code, ik) /= TQS_VALID) return
        geo(1) = g%u
        geo(2:3) = g%z_join
        geo(4) = g%l(3)
        geo(5) = g%a_sq(3)
        geo(6) = g%r_north
        geo(7) = g%r_south
        geo(8) = g%alpha1
        write(*, '(A,A,I0)') name, '_MID_KIND = ', g%mid_kind
        call print_list(name, '_GEOMETRY', geo)

        call tqs_param_neck(p, 5_ik_c, nk(1), nk(2), found, code)
        call assert_int_eq(int(code, ik), TQS_VALID, name // ': neck valid')
        write(*, '(A,A,A)') name, '_NECK_FOUND = ', trim(merge('True ', 'False', found == 1_ik_c))
        call print_list(name, '_NECK', nk)

        call tqs_param_moller_coordinates(p, 5_ik_c, mo(1), mo(2), mo(3), mo(4), mo(5), code)
        call assert_int_eq(int(code, ik), TQS_VALID, name // ': moller valid')
        call print_list(name, '_MOLLER', mo)

        call tqs_param_scission_sigma1(p, 5_ik_c, sc(1), sc(2), defined, code)
        call assert_int_eq(int(code, ik), TQS_VALID, name // ': scission valid')
        call print_list(name, '_SCISSION', sc)

        ! Thetas as pre-rounded literals (pi/8, pi/2, 7pi/8), NOT arithmetic:
        ! Release builds may reassociate 7*PI_C/8 during constant folding;
        ! literals parse identically here and in Python.
        eval_thetas = [0.39269908169872414_rk_c, 1.5707963267948966_rk_c, &
                2.748893571891069_rk_c]
        call tqs_param_radius_and_derivative(p, 5_ik_c, eval_thetas, 3_ik_c, eval_r, eval_dr, code)
        call assert_int_eq(int(code, ik), TQS_VALID, name // ': radius_and_derivative valid')
        call print_list(name, '_DERIV_R', eval_r)
        call print_list(name, '_DERIV_DR', eval_dr)

        call tqs_param_rho_z_grid(p, 5_ik_c, int(N_Z, ik_c), z, rho, drho, code)
        call assert_int_eq(int(code, ik), TQS_VALID, name // ': rho_z_grid valid')
        picked = z(IDX)
        call print_list(name, '_Z', picked)
        picked = rho(IDX)
        call print_list(name, '_RHO', picked)
    end subroutine capture

end program golden_capture
