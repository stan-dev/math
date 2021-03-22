! ------------------------------------------------------------------
! Programmer(s): Cody J. Balos @ LLNL
!                Daniel R. Reynolds @ SMU
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
! This is an example custom NVECTOR based on 2D Fortran arrays.
! ------------------------------------------------------------------

module fnvector_fortran_mod

  use, intrinsic :: iso_c_binding
  use fsundials_nvector_mod

  implicit none

  ! ----------------------------------------------------------------
  type, public :: FVec
    logical                 :: own_data
    integer(c_long)         :: length1
    integer(c_long)         :: length2
    real(c_double), pointer :: data(:,:)
  end type FVec
  ! ----------------------------------------------------------------

contains

  ! ----------------------------------------------------------------
  function FN_VNew_Fortran(n1, n2) result(sunvec_y)

    implicit none
    integer(c_long),    value   :: n1
    integer(c_long),    value   :: n2
    type(N_Vector),     pointer :: sunvec_y
    type(N_Vector_Ops), pointer :: ops
    type(FVec),         pointer :: content

    ! allocate output N_Vector structure
    sunvec_y => FN_VNewEmpty()

    ! allocate and fill content structure
    allocate(content)
    allocate(content%data(n1,n2))
    content%own_data = .true.
    content%length1  = n1
    content%length2  = n2

    ! attach the content structure to the output N_Vector
    sunvec_y%content = c_loc(content)

    ! access the N_Vector ops structure, and set internal function pointers
    call c_f_pointer(sunvec_y%ops, ops)
    ops%nvgetvectorid      = c_funloc(FN_VGetVectorID_Fortran)
    ops%nvdestroy          = c_funloc(FN_VDestroy_Fortran)
    ops%nvgetlength        = c_funloc(FN_VGetLength_Fortran)
    ops%nvconst            = c_funloc(FN_VConst_Fortran)
    ops%nvdotprod          = c_funloc(FN_VDotProd_Fortran)
    ops%nvclone            = c_funloc(FN_VClone_Fortran)
    ops%nvspace            = c_funloc(FN_VSpace_Fortran)
    ops%nvlinearsum        = c_funloc(FN_VLinearSum_Fortran)
    ops%nvprod             = c_funloc(FN_VProd_Fortran)
    ops%nvdiv              = c_funloc(FN_VDiv_Fortran)
    ops%nvscale            = c_funloc(FN_VScale_Fortran)
    ops%nvabs              = c_funloc(FN_VAbs_Fortran)
    ops%nvinv              = c_funloc(FN_VInv_Fortran)
    ops%nvaddconst         = c_funloc(FN_VAddConst_Fortran)
    ops%nvmaxnorm          = c_funloc(FN_VMaxNorm_Fortran)
    ops%nvwrmsnorm         = c_funloc(FN_VWRMSNorm_Fortran)
    ops%nvwrmsnormmask     = c_funloc(FN_VWRMSNormMask_Fortran)
    ops%nvmin              = c_funloc(FN_VMin_Fortran)
    ops%nvwl2norm          = c_funloc(FN_VWL2Norm_Fortran)
    ops%nvl1norm           = c_funloc(FN_VL1Norm_Fortran)
    ops%nvcompare          = c_funloc(FN_VCompare_Fortran)
    ops%nvinvtest          = c_funloc(FN_VInvTest_Fortran)
    ops%nvconstrmask       = c_funloc(FN_VConstrMask_Fortran)
    ops%nvminquotient      = c_funloc(FN_VMinQuotient_Fortran)
    ops%nvdotprodlocal     = c_funloc(FN_VDotProd_Fortran)
    ops%nvmaxnormlocal     = c_funloc(FN_VMaxNorm_Fortran)
    ops%nvminlocal         = c_funloc(FN_VMin_Fortran)
    ops%nvl1normlocal      = c_funloc(FN_VL1Norm_Fortran)
    ops%nvinvtestlocal     = c_funloc(FN_VInvTest_Fortran)
    ops%nvconstrmasklocal  = c_funloc(FN_VConstrMask_Fortran)
    ops%nvminquotientlocal = c_funloc(FN_VMinQuotient_Fortran)
    ops%nvwsqrsumlocal     = c_funloc(FN_VWSqrSum_Fortran)
    ops%nvwsqrsummasklocal = c_funloc(FN_VWSqrSumMask_Fortran)

  end function FN_VNew_Fortran

  ! ----------------------------------------------------------------
  function FN_VMake_Fortran(n1, n2, data) result(sunvec_y)

    implicit none
    integer(c_long),    value   :: n1
    integer(c_long),    value   :: n2
    type(N_Vector),     pointer :: sunvec_y
    type(N_Vector_Ops), pointer :: ops
    type(FVec),         pointer :: content
    real(c_double),     target  :: data(:,:)

    ! allocate output N_Vector structure
    sunvec_y => FN_VNewEmpty()

    ! allocate and fill content structure
    allocate(content)
    content%own_data = .false.
    content%length1  = n1
    content%length2  = n2
    content%data    => data

    ! attach the content structure to the output N_Vector
    sunvec_y%content = c_loc(content)

    ! access the N_Vector ops structure, and set internal function pointers
    call c_f_pointer(sunvec_y%ops, ops)
    ops%nvgetvectorid      = c_funloc(FN_VGetVectorID_Fortran)
    ops%nvdestroy          = c_funloc(FN_VDestroy_Fortran)
    ops%nvgetlength        = c_funloc(FN_VGetLength_Fortran)
    ops%nvconst            = c_funloc(FN_VConst_Fortran)
    ops%nvdotprod          = c_funloc(FN_VDotProd_Fortran)
    ops%nvclone            = c_funloc(FN_VClone_Fortran)
    ops%nvspace            = c_funloc(FN_VSpace_Fortran)
    ops%nvlinearsum        = c_funloc(FN_VLinearSum_Fortran)
    ops%nvprod             = c_funloc(FN_VProd_Fortran)
    ops%nvdiv              = c_funloc(FN_VDiv_Fortran)
    ops%nvscale            = c_funloc(FN_VScale_Fortran)
    ops%nvabs              = c_funloc(FN_VAbs_Fortran)
    ops%nvinv              = c_funloc(FN_VInv_Fortran)
    ops%nvaddconst         = c_funloc(FN_VAddConst_Fortran)
    ops%nvmaxnorm          = c_funloc(FN_VMaxNorm_Fortran)
    ops%nvwrmsnorm         = c_funloc(FN_VWRMSNorm_Fortran)
    ops%nvwrmsnormmask     = c_funloc(FN_VWRMSNormMask_Fortran)
    ops%nvmin              = c_funloc(FN_VMin_Fortran)
    ops%nvwl2norm          = c_funloc(FN_VWL2Norm_Fortran)
    ops%nvl1norm           = c_funloc(FN_VL1Norm_Fortran)
    ops%nvcompare          = c_funloc(FN_VCompare_Fortran)
    ops%nvinvtest          = c_funloc(FN_VInvTest_Fortran)
    ops%nvconstrmask       = c_funloc(FN_VConstrMask_Fortran)
    ops%nvminquotient      = c_funloc(FN_VMinQuotient_Fortran)
    ops%nvdotprodlocal     = c_funloc(FN_VDotProd_Fortran)
    ops%nvmaxnormlocal     = c_funloc(FN_VMaxNorm_Fortran)
    ops%nvminlocal         = c_funloc(FN_VMin_Fortran)
    ops%nvl1normlocal      = c_funloc(FN_VL1Norm_Fortran)
    ops%nvinvtestlocal     = c_funloc(FN_VInvTest_Fortran)
    ops%nvconstrmasklocal  = c_funloc(FN_VConstrMask_Fortran)
    ops%nvminquotientlocal = c_funloc(FN_VMinQuotient_Fortran)
    ops%nvwsqrsumlocal     = c_funloc(FN_VWSqrSum_Fortran)
    ops%nvwsqrsummasklocal = c_funloc(FN_VWSqrSumMask_Fortran)

  end function FN_VMake_Fortran

  ! ----------------------------------------------------------------
  function FN_VGetFVec(sunvec_x) result(x)

    implicit none
    type(N_Vector)         :: sunvec_x
    type(FVec),    pointer :: x

    ! extract Fortran matrix structure to output
    call c_f_pointer(sunvec_x%content, x)

    return

  end function FN_VGetFVec

  ! ----------------------------------------------------------------
  integer(N_Vector_ID) function FN_VGetVectorID_Fortran(sunvec_y) &
    result(id) bind(C)

    implicit none
    type(N_Vector) :: sunvec_y

    id = SUNDIALS_NVEC_CUSTOM
    return

  end function FN_VGetVectorID_Fortran

  ! ----------------------------------------------------------------
  subroutine FN_VDestroy_Fortran(sunvec_y) bind(C)

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

  end subroutine FN_VDestroy_Fortran

  ! ----------------------------------------------------------------
  integer(c_long) function FN_VGetLength_Fortran(sunvec_y) &
      bind(C) result(length)

    implicit none
    type(N_Vector)      :: sunvec_y
    type(FVec), pointer :: y

    y => FN_VGetFVec(sunvec_y)
    length = (y%length1)*(y%length2)
    return

  end function FN_VGetLength_Fortran

  ! ----------------------------------------------------------------
  subroutine FN_VConst_Fortran(const, sunvec_y) bind(C)

    implicit none
    type(N_Vector)          :: sunvec_y
    real(c_double), value   :: const
    type(FVec),     pointer :: y

    ! extract Fortran vector structure to work with
    y => FN_VGetFVec(sunvec_y)

    ! set all entries to const (whole array operation)
    y%data = const
    return

  end subroutine FN_VConst_Fortran

  ! ----------------------------------------------------------------
  real(c_double) function FN_VDotProd_Fortran(sunvec_x, sunvec_y) &
    result(a) bind(C)

    implicit none
    type(N_Vector)         :: sunvec_x, sunvec_y
    type(FVec),    pointer :: x, y

    ! extract Fortran vector structures to work with
    x => FN_VGetFVec(sunvec_x)
    y => FN_VGetFVec(sunvec_y)

    ! do the dot product via Fortran intrinsics
    a = sum(x%data * y%data)
    return

  end function FN_VDotProd_Fortran

  ! ----------------------------------------------------------------
  function FN_VClone_Fortran(sunvec_x) result(y_ptr) bind(C)

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
    allocate(y%data(x%length1,x%length2))
    y%own_data = .true.
    y%length1 = x%length1
    y%length2 = x%length2

    ! attach the content structure to the output N_Vector
    sunvec_y%content = c_loc(y)

    ! set the c_ptr output
    y_ptr = c_loc(sunvec_y)
    return

  end function FN_VClone_Fortran

  ! ----------------------------------------------------------------
  subroutine FN_VSpace_Fortran(sunvec_x, lrw, liw) bind(C)

    implicit none
    type(N_Vector)      :: sunvec_x
    integer(c_int64_t)  :: lrw(1)
    integer(c_int64_t)  :: liw(1)
    type(FVec), pointer :: x

    ! extract Fortran vector structure to work with
    x => FN_VGetFVec(sunvec_x)

    ! set output arguments and return
    lrw(1) = (x%length1)*(x%length2)
    liw(1) = 3
    return

  end subroutine FN_VSpace_Fortran

  ! ----------------------------------------------------------------
  subroutine FN_VLinearSum_Fortran(a, sunvec_x, b, sunvec_y, sunvec_z) &
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

  end subroutine FN_VLinearSum_Fortran

  ! ----------------------------------------------------------------
  subroutine FN_VProd_Fortran(sunvec_x, sunvec_y, sunvec_z) bind(C)

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

  end subroutine FN_VProd_Fortran

  ! ----------------------------------------------------------------
  subroutine FN_VDiv_Fortran(sunvec_x, sunvec_y, sunvec_z) bind(C)

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

  end subroutine FN_VDiv_Fortran

  ! ----------------------------------------------------------------
  subroutine FN_VScale_Fortran(c, sunvec_x, sunvec_z) bind(C)

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

  end subroutine FN_VScale_Fortran

  ! ----------------------------------------------------------------
  subroutine FN_VAbs_Fortran(sunvec_x, sunvec_z) bind(C)

    implicit none
    type(N_Vector)      :: sunvec_x
    type(N_Vector)      :: sunvec_z
    type(FVec), pointer :: x, z

    ! extract Fortran vector structures to work with
    x => FN_VGetFVec(sunvec_x)
    z => FN_VGetFVec(sunvec_z)

    ! perform computation (via whole array ops) and return
    z%data = dabs(x%data)
    return

  end subroutine FN_VAbs_Fortran

  ! ----------------------------------------------------------------
  subroutine FN_VInv_Fortran(sunvec_x, sunvec_z) bind(C)

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

  end subroutine FN_VInv_Fortran

  ! ----------------------------------------------------------------
  subroutine FN_VAddConst_Fortran(sunvec_x, b, sunvec_z) bind(C)

    implicit none
    type(N_Vector)          :: sunvec_x
    real(c_double), value   :: b
    type(N_Vector)          :: sunvec_z
    type(FVec),     pointer :: x, z

    ! extract Fortran vector structures to work with
    x => FN_VGetFVec(sunvec_x)
    z => FN_VGetFVec(sunvec_z)

    ! perform computation (via whole array ops) and return
    z%data = x%data + b
    return

  end subroutine FN_VAddConst_Fortran

  ! ----------------------------------------------------------------
  real(c_double) function FN_VMaxNorm_Fortran(sunvec_x) &
       result(maxnorm) bind(C)

    implicit none
    type(N_Vector)      :: sunvec_x
    type(FVec), pointer :: x

    ! extract Fortran vector structure to work with
    x => FN_VGetFVec(sunvec_x)

    ! perform computation (via whole array ops) and return
    maxnorm = maxval(dabs(x%data))
    return

  end function FN_VMaxNorm_Fortran

  ! ----------------------------------------------------------------
  real(c_double) function FN_VWSqrSum_Fortran(sunvec_x, sunvec_w) &
       result(sqrsum) bind(C)

    implicit none
    type(N_Vector)      :: sunvec_x
    type(N_Vector)      :: sunvec_w
    type(FVec), pointer :: x, w

    ! extract Fortran vector structures to work with
    x => FN_VGetFVec(sunvec_x)
    w => FN_VGetFVec(sunvec_w)

    ! perform computation (via whole array ops) and return
    sqrsum = sum(x%data * x%data * w%data * w%data)
    return

  end function FN_VWSqrSum_Fortran

  ! ----------------------------------------------------------------
  real(c_double) function FN_VWSqrSumMask_Fortran(sunvec_x, sunvec_w, sunvec_id) &
       result(sqrsum) bind(C)

    implicit none
    type(N_Vector)      :: sunvec_x
    type(N_Vector)      :: sunvec_w
    type(N_Vector)      :: sunvec_id
    type(FVec), pointer :: x, w, id
    integer(c_long)     :: i, j

    ! extract Fortran vector structures to work with
    x => FN_VGetFVec(sunvec_x)
    w => FN_VGetFVec(sunvec_w)
    id => FN_VGetFVec(sunvec_id)

    ! perform computation and return
    sqrsum = 0.d0
    do j = 1,x%length2
       do i = 1,x%length1
          if (id%data(i,j) > 0.d0) then
             sqrsum = sqrsum + (x%data(i,j) * w%data(i,j))**2
          end if
       end do
    end do
    return

  end function FN_VWSqrSumMask_Fortran

  ! ----------------------------------------------------------------
  real(c_double) function FN_VWRMSNorm_Fortran(sunvec_x, sunvec_w) &
       result(wrmsnorm) bind(C)

    implicit none
    type(N_Vector)      :: sunvec_x
    type(N_Vector)      :: sunvec_w
    type(FVec), pointer :: x

    ! extract Fortran vector structure to work with
    x => FN_VGetFVec(sunvec_x)

    ! postprocess result from FN_VWSqrSum for result
    wrmsnorm = dsqrt( FN_VWSqrSum_Fortran(sunvec_x, sunvec_w) / (x%length1 * x%length2) )
    return

  end function FN_VWRMSNorm_Fortran

  ! ----------------------------------------------------------------
  real(c_double) function FN_VWRMSNormMask_Fortran(sunvec_x, sunvec_w, sunvec_id) &
       result(wrmsnorm) bind(C)

    implicit none
    type(N_Vector)      :: sunvec_x
    type(N_Vector)      :: sunvec_w
    type(N_Vector)      :: sunvec_id
    type(FVec), pointer :: x

    ! extract Fortran vector structure to work with
    x => FN_VGetFVec(sunvec_x)

    ! postprocess result from FN_VWSqrSumMask for result
    wrmsnorm = dsqrt( FN_VWSqrSumMask_Fortran(sunvec_x, sunvec_w, sunvec_id) / (x%length1 * x%length2) )
    return

  end function FN_VWRMSNormMask_Fortran

  ! ----------------------------------------------------------------
  real(c_double) function FN_VMin_Fortran(sunvec_x) &
       result(mnval) bind(C)

    implicit none
    type(N_Vector)      :: sunvec_x
    type(FVec), pointer :: x

    ! extract Fortran vector structure to work with
    x => FN_VGetFVec(sunvec_x)

    ! perform computation (via whole array ops) and return
    mnval = minval(x%data)
    return

  end function FN_VMin_Fortran

  ! ----------------------------------------------------------------
  real(c_double) function FN_VWL2Norm_Fortran(sunvec_x, sunvec_w) &
       result(wl2norm) bind(C)

    implicit none
    type(N_Vector)      :: sunvec_x
    type(N_Vector)      :: sunvec_w
    type(FVec), pointer :: x

    ! extract Fortran vector structure to work with
    x => FN_VGetFVec(sunvec_x)

    ! postprocess result from FN_VWSqrSum for result
    wl2norm = dsqrt(FN_VWSqrSum_Fortran(sunvec_x, sunvec_w))
    return

  end function FN_VWL2Norm_Fortran

  ! ----------------------------------------------------------------
  real(c_double) function FN_VL1Norm_Fortran(sunvec_x) &
       result(l1norm) bind(C)

    implicit none
    type(N_Vector)      :: sunvec_x
    type(FVec), pointer :: x

    ! extract Fortran vector structure to work with
    x => FN_VGetFVec(sunvec_x)

    ! perform computation (via whole array ops) and return
    l1norm = sum(abs(x%data))
    return

  end function FN_VL1Norm_Fortran

  ! ----------------------------------------------------------------
  subroutine FN_VCompare_Fortran(c, sunvec_x, sunvec_z) bind(C)

    implicit none
    real(c_double), value   :: c
    type(N_Vector)          :: sunvec_x
    type(N_Vector)          :: sunvec_z
    type(FVec),     pointer :: x, z
    integer(c_long)         :: i, j

    ! extract Fortran vector structures to work with
    x => FN_VGetFVec(sunvec_x)
    z => FN_VGetFVec(sunvec_z)

    ! perform operation and return
    do j = 1,x%length2
       do i = 1,x%length1
          if (abs(x%data(i,j)) .ge. c) then
             z%data(i,j) = 1.d0
          else
             z%data(i,j) = 0.d0
          end if
       end do
    end do
    return

  end subroutine FN_VCompare_Fortran

  ! ----------------------------------------------------------------
  integer(c_int) function FN_VInvTest_Fortran(sunvec_x, sunvec_z) &
       result(no_zero_found) bind(C)

    implicit none
    type(N_Vector)      :: sunvec_x
    type(N_Vector)      :: sunvec_z
    type(FVec), pointer :: x, z
    integer(c_long)     :: i, j

    ! extract Fortran vector structures to work with
    x => FN_VGetFVec(sunvec_x)
    z => FN_VGetFVec(sunvec_z)

    ! perform operation and return
    no_zero_found = 1
    do j = 1,x%length2
       do i = 1,x%length1
          if (x%data(i,j) == 0.d0) then
             no_zero_found = 0
          else
             z%data(i,j) = 1.d0 / x%data(i,j)
          end if
       end do
    end do
    return

  end function FN_VInvTest_Fortran

  ! ----------------------------------------------------------------
  integer(c_int) function FN_VConstrMask_Fortran(sunvec_c, sunvec_x, sunvec_m) &
       result(all_good) bind(C)

    implicit none
    type(N_Vector)      :: sunvec_c
    type(N_Vector)      :: sunvec_x
    type(N_Vector)      :: sunvec_m
    type(FVec), pointer :: c, x, m
    integer(c_long)     :: i, j
    logical             :: test

    ! extract Fortran vector structures to work with
    c => FN_VGetFVec(sunvec_c)
    x => FN_VGetFVec(sunvec_x)
    m => FN_VGetFVec(sunvec_m)

    ! perform operation and return
    all_good = 1
    do j = 1,x%length2
       do i = 1,x%length1
          m%data(i,j) = 0.d0

          ! continue if no constraints were set for this variable
          if (c%data(i,j) == 0.d0)  cycle

          ! check if a set constraint has been violated
          test = ((dabs(c%data(i,j)) > 1.5d0 .and. x%data(i,j)*c%data(i,j) .le. 0.d0) .or. &
                  (dabs(c%data(i,j)) > 0.5d0 .and. x%data(i,j)*c%data(i,j) <  0.d0))
          if (test) then
             all_good = 0
             m%data(i,j) = 1.d0
          end if
       end do
    end do
    return

  end function FN_VConstrMask_Fortran

  ! ----------------------------------------------------------------
  real(c_double) function FN_VMinQuotient_Fortran(sunvec_n, sunvec_d) &
       result(minq) bind(C)

    implicit none
    type(N_Vector)      :: sunvec_n
    type(N_Vector)      :: sunvec_d
    type(FVec), pointer :: n, d
    integer(c_long)     :: i, j
    logical             :: notEvenOnce

    ! extract Fortran vector structures to work with
    n => FN_VGetFVec(sunvec_n)
    d => FN_VGetFVec(sunvec_d)

    ! initialize results
    notEvenOnce = .true.
    minq = 1.d307

    ! perform operation and return
    do j = 1,n%length2
       do i = 1,n%length1

          ! skip locations with zero-valued denominator
          if (d%data(i,j) == 0.d0)  cycle

          ! store the first quotient value
          if (notEvenOnce) then
             minq = n%data(i,j)/d%data(i,j)
             notEvenOnce = .false.
          else
             minq = min(minq, n%data(i,j)/d%data(i,j))
          end if
       end do
    end do
    return

  end function FN_VMinQuotient_Fortran

end module fnvector_fortran_mod
! ------------------------------------------------------------------
