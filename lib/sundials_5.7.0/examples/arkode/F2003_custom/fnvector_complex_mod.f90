! ------------------------------------------------------------------
! Programmer(s): Daniel R. Reynolds @ SMU
! ------------------------------------------------------------------
! SUNDIALS Copyright Start
! Copyright (c) 2002-2021, Lawrence Livermore National Security
! and Southern Methodist University.
! All rights reserved.
!
! See the top-level LICENSE and NOTICE files for details.
!
! SPDX-License-Identifier: BSD-3-Clause
! SUNDIALS Copyright End
! ------------------------------------------------------------------
! This is an example custom NVECTOR containing complex data.
! ------------------------------------------------------------------

module fnvector_complex_mod

  use, intrinsic :: iso_c_binding
  use fsundials_nvector_mod

  implicit none

  ! ----------------------------------------------------------------
  type, public :: FVec
    logical                            :: own_data
    integer(c_long)                    :: len
    complex(c_double_complex), pointer :: data(:)
  end type FVec
  ! ----------------------------------------------------------------

contains

  ! ----------------------------------------------------------------
  function FN_VNew_Complex(n) result(sunvec_y)

    implicit none
    integer(c_long),    value   :: n
    type(N_Vector),     pointer :: sunvec_y
    type(N_Vector_Ops), pointer :: ops
    type(FVec),         pointer :: content

    ! allocate output N_Vector structure
    sunvec_y => FN_VNewEmpty()

    ! allocate and fill content structure
    allocate(content)
    allocate(content%data(n))
    content%own_data = .true.
    content%len = n

    ! attach the content structure to the output N_Vector
    sunvec_y%content = c_loc(content)

    ! access the N_Vector ops structure, and set internal function pointers
    call c_f_pointer(sunvec_y%ops, ops)
    ops%nvgetvectorid      = c_funloc(FN_VGetVectorID_Complex)
    ops%nvdestroy          = c_funloc(FN_VDestroy_Complex)
    ops%nvgetlength        = c_funloc(FN_VGetLength_Complex)
    ops%nvconst            = c_funloc(FN_VConst_Complex)
    ops%nvclone            = c_funloc(FN_VClone_Complex)
    ops%nvspace            = c_funloc(FN_VSpace_Complex)
    ops%nvlinearsum        = c_funloc(FN_VLinearSum_Complex)
    ops%nvprod             = c_funloc(FN_VProd_Complex)
    ops%nvdiv              = c_funloc(FN_VDiv_Complex)
    ops%nvscale            = c_funloc(FN_VScale_Complex)
    ops%nvabs              = c_funloc(FN_VAbs_Complex)
    ops%nvinv              = c_funloc(FN_VInv_Complex)
    ops%nvaddconst         = c_funloc(FN_VAddConst_Complex)
    ops%nvmaxnorm          = c_funloc(FN_VMaxNorm_Complex)
    ops%nvwrmsnorm         = c_funloc(FN_VWRMSNorm_Complex)
    ops%nvwrmsnormmask     = c_funloc(FN_VWRMSNormMask_Complex)
    ops%nvmin              = c_funloc(FN_VMin_Complex)
    ops%nvwl2norm          = c_funloc(FN_VWL2Norm_Complex)
    ops%nvl1norm           = c_funloc(FN_VL1Norm_Complex)
    ops%nvinvtest          = c_funloc(FN_VInvTest_Complex)
    ops%nvmaxnormlocal     = c_funloc(FN_VMaxNorm_Complex)
    ops%nvminlocal         = c_funloc(FN_VMin_Complex)
    ops%nvl1normlocal      = c_funloc(FN_VL1Norm_Complex)
    ops%nvinvtestlocal     = c_funloc(FN_VInvTest_Complex)
    ops%nvwsqrsumlocal     = c_funloc(FN_VWSqrSum_Complex)
    ops%nvwsqrsummasklocal = c_funloc(FN_VWSqrSumMask_Complex)

  end function FN_VNew_Complex

  ! ----------------------------------------------------------------
  function FN_VMake_Complex(n, data) result(sunvec_y)

    implicit none
    integer(c_long),           value   :: n
    type(N_Vector),            pointer :: sunvec_y
    type(N_Vector_Ops),        pointer :: ops
    type(FVec),                pointer :: content
    complex(c_double_complex), target  :: data(:)

    ! allocate output N_Vector structure
    sunvec_y => FN_VNewEmpty()

    ! allocate and fill content structure
    allocate(content)
    content%own_data = .false.
    content%len = n
    content%data => data

    ! attach the content structure to the output N_Vector
    sunvec_y%content = c_loc(content)

    ! access the N_Vector ops structure, and set internal function pointers
    call c_f_pointer(sunvec_y%ops, ops)
    ops%nvgetvectorid      = c_funloc(FN_VGetVectorID_Complex)
    ops%nvdestroy          = c_funloc(FN_VDestroy_Complex)
    ops%nvgetlength        = c_funloc(FN_VGetLength_Complex)
    ops%nvconst            = c_funloc(FN_VConst_Complex)
    ops%nvclone            = c_funloc(FN_VClone_Complex)
    ops%nvspace            = c_funloc(FN_VSpace_Complex)
    ops%nvlinearsum        = c_funloc(FN_VLinearSum_Complex)
    ops%nvprod             = c_funloc(FN_VProd_Complex)
    ops%nvdiv              = c_funloc(FN_VDiv_Complex)
    ops%nvscale            = c_funloc(FN_VScale_Complex)
    ops%nvabs              = c_funloc(FN_VAbs_Complex)
    ops%nvinv              = c_funloc(FN_VInv_Complex)
    ops%nvaddconst         = c_funloc(FN_VAddConst_Complex)
    ops%nvmaxnorm          = c_funloc(FN_VMaxNorm_Complex)
    ops%nvwrmsnorm         = c_funloc(FN_VWRMSNorm_Complex)
    ops%nvwrmsnormmask     = c_funloc(FN_VWRMSNormMask_Complex)
    ops%nvmin              = c_funloc(FN_VMin_Complex)
    ops%nvwl2norm          = c_funloc(FN_VWL2Norm_Complex)
    ops%nvl1norm           = c_funloc(FN_VL1Norm_Complex)
    ops%nvinvtest          = c_funloc(FN_VInvTest_Complex)
    ops%nvmaxnormlocal     = c_funloc(FN_VMaxNorm_Complex)
    ops%nvminlocal         = c_funloc(FN_VMin_Complex)
    ops%nvl1normlocal      = c_funloc(FN_VL1Norm_Complex)
    ops%nvinvtestlocal     = c_funloc(FN_VInvTest_Complex)
    ops%nvwsqrsumlocal     = c_funloc(FN_VWSqrSum_Complex)
    ops%nvwsqrsummasklocal = c_funloc(FN_VWSqrSumMask_Complex)

  end function FN_VMake_Complex

  ! ----------------------------------------------------------------
  function FN_VGetFVec(sunvec_x) result(x)

    implicit none
    type(N_Vector) :: sunvec_x
    type(FVec), pointer :: x

    ! extract Fortran matrix structure to output
    call c_f_pointer(sunvec_x%content, x)

    return

  end function FN_VGetFVec

  ! ----------------------------------------------------------------
  integer(N_Vector_ID) function FN_VGetVectorID_Complex(sunvec_y) &
    result(id) bind(C)

    implicit none
    type(N_Vector) :: sunvec_y

    id = SUNDIALS_NVEC_CUSTOM
    return

  end function FN_VGetVectorID_Complex

  ! ----------------------------------------------------------------
  subroutine FN_VDestroy_Complex(sunvec_y) bind(C)

    implicit none
    type(N_Vector), target  :: sunvec_y
    type(FVec),     pointer :: y

    ! access FVec structure
    y => FN_VGetFVec(sunvec_y)

    ! if vector owns the data, then deallocate
    if (y%own_data)  deallocate(y%data)

    ! deallocate the underlying Fortran object (the content)
    deallocate(y)

    ! set N_Vector structure members to NULL and return
    sunvec_y%content = C_NULL_PTR

    ! deallocate overall N_Vector structure
    call FN_VFreeEmpty(sunvec_y)

    return

  end subroutine FN_VDestroy_Complex

  ! ----------------------------------------------------------------
  integer(c_long) function FN_VGetLength_Complex(sunvec_y) &
      bind(C) result(length)

    implicit none
    type(N_Vector)         :: sunvec_y
    type(FVec),    pointer :: y

    y => FN_VGetFVec(sunvec_y)
    length = y%len
    return

  end function FN_VGetLength_Complex

  ! ----------------------------------------------------------------
  subroutine FN_VConst_Complex(const, sunvec_y) bind(C)

    implicit none
    type(N_Vector)          :: sunvec_y
    real(c_double), value   :: const
    type(FVec),     pointer :: y

    ! extract Fortran vector structure to work with
    y => FN_VGetFVec(sunvec_y)

    ! set all entries to const (whole array operation)
    y%data = dcmplx(const, 0.d0)
    return

  end subroutine FN_VConst_Complex

  ! ----------------------------------------------------------------
  function FN_VClone_Complex(sunvec_x) result(y_ptr) bind(C)

    implicit none
    type(N_Vector)          :: sunvec_x
    type(N_Vector), pointer :: sunvec_y
    type(c_ptr)             :: y_ptr
    integer(c_int)          :: retval
    type(FVec),     pointer :: x, y

    ! extract Fortran vector structure to work with
    x => FN_VGetFVec(sunvec_x)

    ! allocate output N_Vector structure
    sunvec_y => FN_VNewEmpty()

    ! copy operations from x into y
    retval = FN_VCopyOps(sunvec_x, sunvec_y)

    ! allocate and clone content structure
    allocate(y)
    allocate(y%data(x%len))
    y%own_data = .true.
    y%len = x%len

    ! attach the content structure to the output N_Vector
    sunvec_y%content = c_loc(y)

    ! set the c_ptr output
    y_ptr = c_loc(sunvec_y)
    return

  end function FN_VClone_Complex

  ! ----------------------------------------------------------------
  subroutine FN_VSpace_Complex(sunvec_x, lrw, liw) bind(C)

    implicit none
    type(N_Vector)      :: sunvec_x
    integer(c_int64_t)  :: lrw(1)
    integer(c_int64_t)  :: liw(1)
    type(FVec), pointer :: x

    ! extract Fortran vector structure to work with
    x => FN_VGetFVec(sunvec_x)

    ! set output arguments and return (multiply lrw by 2 since complex)
    lrw(1) = 2*x%len
    liw(1) = 3
    return

  end subroutine FN_VSpace_Complex

  ! ----------------------------------------------------------------
  subroutine FN_VLinearSum_Complex(a, sunvec_x, b, sunvec_y, sunvec_z) &
       bind(C)

    implicit none
    type(N_Vector)          :: sunvec_x
    type(N_Vector)          :: sunvec_y
    type(N_Vector)          :: sunvec_z
    real(c_double), value   :: a
    real(c_double), value   :: b
    type(FVec),     pointer :: x, y, z

    ! extract Fortran vector structures to work with
    x => FN_VGetFVec(sunvec_x)
    y => FN_VGetFVec(sunvec_y)
    z => FN_VGetFVec(sunvec_z)

    ! perform computation (via whole array ops) and return
    z%data = a*x%data + b*y%data
    return

  end subroutine FN_VLinearSum_Complex

  ! ----------------------------------------------------------------
  subroutine FN_VProd_Complex(sunvec_x, sunvec_y, sunvec_z) bind(C)

    implicit none
    type(N_Vector)      :: sunvec_x
    type(N_Vector)      :: sunvec_y
    type(N_Vector)      :: sunvec_z
    type(FVec), pointer :: x, y, z

    ! extract Fortran vector structures to work with
    x => FN_VGetFVec(sunvec_x)
    y => FN_VGetFVec(sunvec_y)
    z => FN_VGetFVec(sunvec_z)

    ! perform computation (via whole array ops) and return
    z%data = x%data * y%data
    return

  end subroutine FN_VProd_Complex

  ! ----------------------------------------------------------------
  subroutine FN_VDiv_Complex(sunvec_x, sunvec_y, sunvec_z) bind(C)

    implicit none
    type(N_Vector)      :: sunvec_x
    type(N_Vector)      :: sunvec_y
    type(N_Vector)      :: sunvec_z
    type(FVec), pointer :: x, y, z

    ! extract Fortran vector structures to work with
    x => FN_VGetFVec(sunvec_x)
    y => FN_VGetFVec(sunvec_y)
    z => FN_VGetFVec(sunvec_z)

    ! perform computation (via whole array ops) and return
    z%data = x%data / y%data
    return

  end subroutine FN_VDiv_Complex

  ! ----------------------------------------------------------------
  subroutine FN_VScale_Complex(c, sunvec_x, sunvec_z) bind(C)

    implicit none
    real(c_double), value   :: c
    type(N_Vector)          :: sunvec_x
    type(N_Vector)          :: sunvec_z
    type(FVec),     pointer :: x, z

    ! extract Fortran vector structures to work with
    x => FN_VGetFVec(sunvec_x)
    z => FN_VGetFVec(sunvec_z)

    ! perform computation (via whole array ops) and return
    z%data = c * x%data
    return

  end subroutine FN_VScale_Complex

  ! ----------------------------------------------------------------
  subroutine FN_VAbs_Complex(sunvec_x, sunvec_z) bind(C)

    implicit none
    type(N_Vector)      :: sunvec_x
    type(N_Vector)      :: sunvec_z
    type(FVec), pointer :: x, z

    ! extract Fortran vector structures to work with
    x => FN_VGetFVec(sunvec_x)
    z => FN_VGetFVec(sunvec_z)

    ! perform computation (via whole array ops) and return
    z%data = abs(x%data)
    return

  end subroutine FN_VAbs_Complex

  ! ----------------------------------------------------------------
  subroutine FN_VInv_Complex(sunvec_x, sunvec_z) bind(C)

    implicit none
    type(N_Vector)      :: sunvec_x
    type(N_Vector)      :: sunvec_z
    type(FVec), pointer :: x, z

    ! extract Fortran vector structures to work with
    x => FN_VGetFVec(sunvec_x)
    z => FN_VGetFVec(sunvec_z)

    ! perform computation (via whole array ops) and return
    z%data = 1.d0 / x%data
    return

  end subroutine FN_VInv_Complex

  ! ----------------------------------------------------------------
  subroutine FN_VAddConst_Complex(sunvec_x, b, sunvec_z) bind(C)

    implicit none
    type(N_Vector)          :: sunvec_x
    real(c_double), value   :: b
    type(N_Vector)          :: sunvec_z
    type(FVec),     pointer :: x, z

    ! extract Fortran vector structures to work with
    x => FN_VGetFVec(sunvec_x)
    z => FN_VGetFVec(sunvec_z)

    ! perform computation (via whole array ops) and return
    z%data = x%data + dcmplx(b, 0.d0)
    return

  end subroutine FN_VAddConst_Complex

  ! ----------------------------------------------------------------
  real(c_double) function FN_VMaxNorm_Complex(sunvec_x) &
       result(maxnorm) bind(C)

    implicit none
    type(N_Vector)      :: sunvec_x
    type(FVec), pointer :: x

    ! extract Fortran vector structure to work with
    x => FN_VGetFVec(sunvec_x)

    ! perform computation (via whole array ops) and return
    maxnorm = maxval(abs(x%data))
    return

  end function FN_VMaxNorm_Complex

  ! ----------------------------------------------------------------
  real(c_double) function FN_VWSqrSum_Complex(sunvec_x, sunvec_w) &
       result(sqrsum) bind(C)

    implicit none
    type(N_Vector)      :: sunvec_x
    type(N_Vector)      :: sunvec_w
    type(FVec), pointer :: x, w

    ! extract Fortran vector structures to work with
    x => FN_VGetFVec(sunvec_x)
    w => FN_VGetFVec(sunvec_w)

    ! perform computation (via whole array ops) and return
    sqrsum = sum(abs(x%data)**2 * abs(w%data)**2)
    return

  end function FN_VWSqrSum_Complex

  ! ----------------------------------------------------------------
  real(c_double) function FN_VWSqrSumMask_Complex(sunvec_x, sunvec_w, sunvec_id) &
       result(sqrsum) bind(C)

    implicit none
    type(N_Vector)      :: sunvec_x
    type(N_Vector)      :: sunvec_w
    type(N_Vector)      :: sunvec_id
    type(FVec), pointer :: x, w, id
    integer(c_long)     :: i

    ! extract Fortran vector structures to work with
    x => FN_VGetFVec(sunvec_x)
    w => FN_VGetFVec(sunvec_w)
    id => FN_VGetFVec(sunvec_id)

    ! perform computation and return
    sqrsum = 0.d0
    do i = 1,x%len
       if (real(id%data(i)) > 0.d0) then
          sqrsum = sqrsum + (abs(x%data(i)) * abs(w%data(i)))**2
       end if
    end do
    return

  end function FN_VWSqrSumMask_Complex

  ! ----------------------------------------------------------------
  real(c_double) function FN_VWRMSNorm_Complex(sunvec_x, sunvec_w) &
       result(wrmsnorm) bind(C)

    implicit none
    type(N_Vector)      :: sunvec_x
    type(N_Vector)      :: sunvec_w
    type(FVec), pointer :: x

    ! extract Fortran vector structure to work with
    x => FN_VGetFVec(sunvec_x)

    ! postprocess result from FN_VWSqrSum for result
    wrmsnorm = dsqrt( FN_VWSqrSum_Complex(sunvec_x, sunvec_w) / x%len )
    return

  end function FN_VWRMSNorm_Complex

  ! ----------------------------------------------------------------
  real(c_double) function FN_VWRMSNormMask_Complex(sunvec_x, sunvec_w, sunvec_id) &
       result(wrmsnorm) bind(C)

    implicit none
    type(N_Vector)      :: sunvec_x
    type(N_Vector)      :: sunvec_w
    type(N_Vector)      :: sunvec_id
    type(FVec), pointer :: x

    ! extract Fortran vector structure to work with
    x => FN_VGetFVec(sunvec_x)

    ! postprocess result from FN_VWSqrSumMask for result
    wrmsnorm = dsqrt( FN_VWSqrSumMask_Complex(sunvec_x, sunvec_w, sunvec_id) / x%len )
    return

  end function FN_VWRMSNormMask_Complex

  ! ----------------------------------------------------------------
  real(c_double) function FN_VMin_Complex(sunvec_x) &
       result(mnval) bind(C)

    implicit none
    type(N_Vector)      :: sunvec_x
    type(FVec), pointer :: x

    ! extract Fortran vector structure to work with
    x => FN_VGetFVec(sunvec_x)

    ! perform computation (via whole array ops) and return
    ! Note: SUNDIALS only wants the minimum of all real components
    mnval = minval(real(x%data))
    return

  end function FN_VMin_Complex

  ! ----------------------------------------------------------------
  real(c_double) function FN_VWL2Norm_Complex(sunvec_x, sunvec_w) &
       result(wl2norm) bind(C)

    implicit none
    type(N_Vector)      :: sunvec_x
    type(N_Vector)      :: sunvec_w
    type(FVec), pointer :: x

    ! extract Fortran vector structure to work with
    x => FN_VGetFVec(sunvec_x)

    ! postprocess result from FN_VWSqrSum for result
    wl2norm = dsqrt(FN_VWSqrSum_Complex(sunvec_x, sunvec_w))
    return

  end function FN_VWL2Norm_Complex

  ! ----------------------------------------------------------------
  real(c_double) function FN_VL1Norm_Complex(sunvec_x) &
       result(l1norm) bind(C)

    implicit none
    type(N_Vector)      :: sunvec_x
    type(FVec), pointer :: x

    ! extract Fortran vector structure to work with
    x => FN_VGetFVec(sunvec_x)

    ! perform computation (via whole array ops) and return
    l1norm = sum(abs(x%data))
    return

  end function FN_VL1Norm_Complex

  ! ----------------------------------------------------------------
  integer(c_int) function FN_VInvTest_Complex(sunvec_x, sunvec_z) &
       result(no_zero_found) bind(C)

    implicit none
    type(N_Vector)      :: sunvec_x
    type(N_Vector)      :: sunvec_z
    type(FVec), pointer :: x, z
    integer(c_long)     :: i

    ! extract Fortran vector structures to work with
    x => FN_VGetFVec(sunvec_x)
    z => FN_VGetFVec(sunvec_z)

    ! perform operation and return
    no_zero_found = 1
    do i = 1,x%len
       if (x%data(i) == dcmplx(0.d0, 0.d0)) then
          no_zero_found = 0
       else
          z%data(i) = 1.d0 / x%data(i)
       end if
    end do
    return

  end function FN_VInvTest_Complex

end module fnvector_complex_mod
! ------------------------------------------------------------------
