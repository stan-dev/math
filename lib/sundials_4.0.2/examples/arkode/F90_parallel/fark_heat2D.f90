!-----------------------------------------------------------------
! Programmer(s): Daniel R. Reynolds @ SMU
!-----------------------------------------------------------------
! SUNDIALS Copyright Start
! Copyright (c) 2002-2019, Lawrence Livermore National Security
! and Southern Methodist University.
! All rights reserved.
!
! See the top-level LICENSE and NOTICE files for details.
!
! SPDX-License-Identifier: BSD-3-Clause
! SUNDIALS Copyright End
!-----------------------------------------------------------------
! Example problem:
! 
! The following test simulates a simple anisotropic 2D heat 
! equation,
!    u_t = kx*u_xx + ky*u_yy + h,
! for t in [0, 10], (x,y) in [0, 1]^2, with initial conditions
!    u(0,x,y) =  0,
! stationary boundary conditions, i.e. 
!    u_t(t,0,y) = u_t(t,1,y) = u_t(t,x,0) = u_t(t,x,1) = 0,
! and a heat source of the form
!    h(x,y) = sin(pi*x)*sin(2*pi*y).
!
! Under this setup, the problem has an analytical solution:
!    u(t,x,y) = a(t)*sin(pi*x)*sin(2*pi*y), where
!    a(t) = (1 - exp(-(kx+4*ky)*pi^2*t)) / ((kx+4*ky)*pi^2).
! 
! The spatial derivatives are computed using second-order 
! centered differences, with the data distributed over nx*ny
! points on a uniform spatial grid.
!
! The spatial grid parameters nx and ny, the parameters kx and ky, 
! as well as the desired relative and absolute solver tolerances, 
! are provided in the input file input_heat2D.txt.
! 
! This program solves the problem with a DIRK method.  This 
! employs a Newton iteration with the PCG iterative linear solver, 
! which itself uses a Jacobi preconditioner.  The example uses the 
! built-in finite-difference Jacobian-vector product routine, but 
! supplies both the RHS and preconditioner setup/solve functions.
!
! 20 outputs are printed at equal intervals, and run statistics 
! are printed at the end.
!-----------------------------------------------------------------

!-----------------------------------------------------------------
! Main driver program
!-----------------------------------------------------------------
module UserData
  !---------------------------------------------------------------
  ! Description: 
  !    Module containing problem-defining parameters, as well as
  !    data buffers for MPI exchanges with neighboring processes.
  !    Also contains routines to:
  !      (a) initialize the module
  !      (b) perform exchanges
  !      (c) free module data.
  !---------------------------------------------------------------
  implicit none
  include "sundials/sundials_fconfig.h"
  save

  integer*8 :: nx          ! global number of x grid points 
  integer*8 :: ny          ! global number of y grid points 
  integer*8 :: is          ! global x indices of this subdomain
  integer*8 :: ie
  integer*8 :: js          ! global y indices of this subdomain
  integer*8 :: je
  integer*8 :: nxl         ! local number of x grid points 
  integer*8 :: nyl         ! local number of y grid points 
  real*8    :: dx          ! x-directional mesh spacing 
  real*8    :: dy          ! y-directional mesh spacing 
  real*8    :: kx          ! x-directional diffusion coefficient 
  real*8    :: ky          ! y-directional diffusion coefficient 
  real(kind=REALTYPE), dimension(:,:), allocatable :: h    ! heat source vector
  real(kind=REALTYPE), dimension(:,:), allocatable :: d    ! inverse of Jacobian diagonal
  integer :: comm                             ! communicator object
  integer :: myid                             ! MPI process ID
  integer :: nprocs                           ! total number of MPI processes
  logical :: HaveNbor(2,2)                    ! flags denoting neighbor on boundary
  real(kind=REALTYPE), dimension(:), allocatable :: Erecv  ! receive buffers for neighbor exchange
  real(kind=REALTYPE), dimension(:), allocatable :: Wrecv
  real(kind=REALTYPE), dimension(:), allocatable :: Nrecv
  real(kind=REALTYPE), dimension(:), allocatable :: Srecv
  real(kind=REALTYPE), dimension(:), allocatable :: Esend  ! send buffers for neighbor exchange
  real(kind=REALTYPE), dimension(:), allocatable :: Wsend
  real(kind=REALTYPE), dimension(:), allocatable :: Nsend
  real(kind=REALTYPE), dimension(:), allocatable :: Ssend

contains

  !---------------------------------------------------------------
  ! Initialize memory allocated within Userdata (set to defaults)
  !---------------------------------------------------------------
  subroutine InitUserData()
    implicit none
    include "mpif.h"
    nx = 0
    ny = 0
    is = 0
    ie = 0
    js = 0
    je = 0
    nxl = 0
    nyl = 0
    dx = 0.d0
    dy = 0.d0
    kx = 0.d0
    ky = 0.d0
    if (allocated(h))  deallocate(h)
    if (allocated(d))  deallocate(d)
    comm = MPI_COMM_WORLD
    myid = 0
    nprocs = 0
    HaveNbor = .false.
    if (allocated(Erecv))  deallocate(Erecv)
    if (allocated(Wrecv))  deallocate(Wrecv)
    if (allocated(Nrecv))  deallocate(Nrecv)
    if (allocated(Srecv))  deallocate(Srecv)
    if (allocated(Esend))  deallocate(Esend)
    if (allocated(Wsend))  deallocate(Wsend)
    if (allocated(Nsend))  deallocate(Nsend)
    if (allocated(Ssend))  deallocate(Ssend)
  end subroutine InitUserData
  !---------------------------------------------------------------


  !---------------------------------------------------------------
  ! Set up parallel decomposition
  !---------------------------------------------------------------
  subroutine SetupDecomp(ierr)
    ! declarations
    implicit none
    include "mpif.h"
    integer, intent(out) :: ierr
    integer :: dims(2), periods(2), coords(2)
    
    ! internals

    ! check that this has not been called before
    if (allocated(h) .or. allocated(d)) then
       write(0,*) "SetupDecomp warning: parallel decomposition already set up"
       ierr = 1
       return
    end if

    ! get suggested parallel decomposition
    dims = (/0, 0/)
    call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)
    if (ierr /= MPI_SUCCESS) then
       write(0,*) "Error in MPI_Comm_size = " , ierr
       return
    end if
    call MPI_Dims_create(nprocs, 2, dims, ierr)
    if (ierr /= MPI_SUCCESS) then
       write(0,*) "Error in MPI_Dims_create = " , ierr
       return
    end if

    ! set up 2D Cartesian communicator
    periods = (/0, 0/)
    call MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, comm, ierr)
    if (ierr /= MPI_SUCCESS) then
       write(0,*) "Error in MPI_Cart_create = " , ierr
       return
    end if
    call MPI_Comm_rank(comm, myid, ierr)
    if (ierr /= MPI_SUCCESS) then
       write(0,*) "Error in MPI_Comm_rank = " , ierr
       return
    end if

    ! determine local extents
    call MPI_Cart_get(comm, 2, dims, periods, coords, ierr)
    if (ierr /= MPI_SUCCESS) then
       write(0,*) "Error in MPI_Cart_get = " , ierr
       return
    end if
    is = nx*coords(1)/dims(1) + 1
    ie = nx*(coords(1)+1)/dims(1)
    js = ny*coords(2)/dims(2) + 1
    je = ny*(coords(2)+1)/dims(2)
    nxl = ie-is+1
    nyl = je-js+1

    ! determine if I have neighbors, and allocate exchange buffers
    HaveNbor(1,1) = (is /= 1)
    HaveNbor(1,2) = (ie /= nx)
    HaveNbor(2,1) = (js /= 1)
    HaveNbor(2,2) = (je /= ny)
    if (HaveNbor(1,1)) then
       allocate(Wrecv(nyl))
       allocate(Wsend(nyl))
    endif
    if (HaveNbor(1,2)) then
       allocate(Erecv(nyl))
       allocate(Esend(nyl))
    endif
    if (HaveNbor(2,1)) then
       allocate(Srecv(nxl))
       allocate(Ssend(nxl))
    endif
    if (HaveNbor(2,2)) then
       allocate(Nrecv(nxl))
       allocate(Nsend(nxl))
    endif

    ! allocate temporary vectors
    allocate(h(nxl,nyl))    ! Create vector for heat source
    allocate(d(nxl,nyl))    ! Create vector for Jacobian diagonal
    
    ierr = 0     ! return with success flag
    return
  end subroutine SetupDecomp
  !---------------------------------------------------------------


  !---------------------------------------------------------------
  ! Perform neighbor exchange
  !---------------------------------------------------------------
  subroutine Exchange(y, ierr)
    ! declarations
    implicit none
    include "mpif.h"
    real(kind=REALTYPE),  intent(in)  :: y(nxl,nyl)
    integer, intent(out) :: ierr
    integer :: reqSW, reqSE, reqSS, reqSN, reqRW, reqRE, reqRS, reqRN;
    integer :: stat(MPI_STATUS_SIZE)
    integer*8 :: i
    integer :: ipW, ipE, ipS, ipN
    integer :: coords(2), dims(2), periods(2), nbcoords(2)
    
    ! internals

    ! MPI neighborhood information
    call MPI_Cart_get(comm, 2, dims, periods, coords, ierr)
    if (ierr /= MPI_SUCCESS) then
       write(0,*) "Error in MPI_Cart_get = ", ierr
       return
    endif
    if (HaveNbor(1,1)) then
       nbcoords = (/ coords(1)-1, coords(2) /)
       call MPI_Cart_rank(comm, nbcoords, ipW, ierr)
       if (ierr /= MPI_SUCCESS) then 
          write(0,*) "Error in MPI_Cart_rank = ", ierr
          return
       endif
    endif
    if (HaveNbor(1,2)) then
       nbcoords = (/ coords(1)+1, coords(2) /)
       call MPI_Cart_rank(comm, nbcoords, ipE, ierr)
       if (ierr /= MPI_SUCCESS) then 
          write(0,*) "Error in MPI_Cart_rank = ", ierr
          return
       endif
    endif
    if (HaveNbor(2,1)) then
       nbcoords = (/ coords(1), coords(2)-1 /)
       call MPI_Cart_rank(comm, nbcoords, ipS, ierr)
       if (ierr /= MPI_SUCCESS) then 
          write(0,*) "Error in MPI_Cart_rank = ", ierr
          return
       endif
    endif
    if (HaveNbor(2,2)) then
       nbcoords = (/ coords(1), coords(2)+1 /)
       call MPI_Cart_rank(comm, nbcoords, ipN, ierr)
       if (ierr /= MPI_SUCCESS) then 
          write(0,*) "Error in MPI_Cart_rank = ", ierr
          return
       endif
    endif

    ! open Irecv buffers
    if (HaveNbor(1,1)) then
       call MPI_Irecv(Wrecv, nyl, MPI_DOUBLE_PRECISION, ipW, &
                      MPI_ANY_TAG, comm, reqRW, ierr)
       if (ierr /= MPI_SUCCESS) then 
          write(0,*) "Error in MPI_Irecv = ", ierr
          return
       endif
    endif
    if (HaveNbor(1,2)) then
       call MPI_Irecv(Erecv, nyl, MPI_DOUBLE_PRECISION, ipE, &
                      MPI_ANY_TAG, comm, reqRE, ierr)
       if (ierr /= MPI_SUCCESS) then 
          write(0,*) "Error in MPI_Irecv = ", ierr
          return
       endif
    endif
    if (HaveNbor(2,1)) then
       call MPI_Irecv(Srecv, nxl, MPI_DOUBLE_PRECISION, ipS, &
                      MPI_ANY_TAG, comm, reqRS, ierr)
       if (ierr /= MPI_SUCCESS) then 
          write(0,*) "Error in MPI_Irecv = ", ierr
          return
       endif
    endif
    if (HaveNbor(2,2)) then
       call MPI_Irecv(Nrecv, nxl, MPI_DOUBLE_PRECISION, ipN, &
                      MPI_ANY_TAG, comm, reqRN, ierr)
       if (ierr /= MPI_SUCCESS) then 
          write(0,*) "Error in MPI_Irecv = ", ierr
          return
       endif
    endif

    ! send data
    if (HaveNbor(1,1)) then
       do i=1,nyl
          Wsend(i) = y(1,i)
       enddo
       call MPI_Isend(Wsend, nyl, MPI_DOUBLE_PRECISION, ipW, 0, &
                      comm, reqSW, ierr)
       if (ierr /= MPI_SUCCESS) then 
          write(0,*) "Error in MPI_Isend = ", ierr
          return
       endif
    endif
    if (HaveNbor(1,2)) then
       do i=1,nyl
          Esend(i) = y(nxl,i)
       enddo
       call MPI_Isend(Esend, nyl, MPI_DOUBLE_PRECISION, ipE, 1, &
                      comm, reqSE, ierr)
       if (ierr /= MPI_SUCCESS) then 
          write(0,*) "Error in MPI_Isend = ", ierr
          return
       endif
    endif
    if (HaveNbor(2,1)) then
       do i=1,nxl
          Ssend(i) = y(i,1)
       enddo
       call MPI_Isend(Ssend, nxl, MPI_DOUBLE_PRECISION, ipS, 2, &
                      comm, reqSS, ierr)
       if (ierr /= MPI_SUCCESS) then 
          write(0,*) "Error in MPI_Isend = ", ierr
          return
       endif
    endif
    if (HaveNbor(2,2)) then
       do i=1,nxl
          Nsend(i) = y(i,nyl)
       enddo
       call MPI_Isend(Nsend, nxl, MPI_DOUBLE_PRECISION, ipN, 3, &
                      comm, reqSN, ierr)
       if (ierr /= MPI_SUCCESS) then 
          write(0,*) "Error in MPI_Isend = ", ierr
          return
       endif
    endif

    ! wait for messages to finish
    if (HaveNbor(1,1)) then
       call MPI_Wait(reqRW, stat, ierr)
       if (ierr /= MPI_SUCCESS) then 
          write(0,*) "Error in MPI_Wait = ", ierr
          return
       endif
       call MPI_Wait(reqSW, stat, ierr)
       if (ierr /= MPI_SUCCESS) then 
          write(0,*) "Error in MPI_Wait = ", ierr
          return
       endif
    endif
    if (HaveNbor(1,2)) then
       call MPI_Wait(reqRE, stat, ierr)
       if (ierr /= MPI_SUCCESS) then 
          write(0,*) "Error in MPI_Wait = ", ierr
          return
       endif
       call MPI_Wait(reqSE, stat, ierr)
       if (ierr /= MPI_SUCCESS) then 
          write(0,*) "Error in MPI_Wait = ", ierr
          return
       endif
    endif
    if (HaveNbor(2,1)) then
       call MPI_Wait(reqRS, stat, ierr)
       if (ierr /= MPI_SUCCESS) then 
          write(0,*) "Error in MPI_Wait = ", ierr
          return
       endif
       call MPI_Wait(reqSS, stat, ierr)
       if (ierr /= MPI_SUCCESS) then 
          write(0,*) "Error in MPI_Wait = ", ierr
          return
       endif
    endif
    if (HaveNbor(2,2)) then
       call MPI_Wait(reqRN, stat, ierr)
       if (ierr /= MPI_SUCCESS) then 
          write(0,*) "Error in MPI_Wait = ", ierr
          return
       endif
       call MPI_Wait(reqSN, stat, ierr)
       if (ierr /= MPI_SUCCESS) then 
          write(0,*) "Error in MPI_Wait = ", ierr
          return
       endif
    endif

    ierr = 0     ! return with success flag
    return
  end subroutine Exchange
  !---------------------------------------------------------------


  !---------------------------------------------------------------
  ! Free memory allocated within Userdata
  !---------------------------------------------------------------
  subroutine FreeUserData(ierr)
    implicit none
    integer, intent(out) :: ierr
    if (allocated(h))      deallocate(h)
    if (allocated(d))      deallocate(d)
    if (allocated(Wrecv))  deallocate(Wrecv)
    if (allocated(Wsend))  deallocate(Wsend)
    if (allocated(Erecv))  deallocate(Erecv)
    if (allocated(Esend))  deallocate(Esend)
    if (allocated(Srecv))  deallocate(Srecv)
    if (allocated(Ssend))  deallocate(Ssend)
    if (allocated(Nrecv))  deallocate(Nrecv)
    if (allocated(Nsend))  deallocate(Nsend)
    ierr = 0     ! return with success flag
    return
  end subroutine FreeUserData
  !---------------------------------------------------------------


  !---------------------------------------------------------------
  ! RMS norm function for parallel array data
  !---------------------------------------------------------------
  subroutine PRMS(y,yrms,ierr)
    ! declarations
    implicit none
    include "mpif.h"
    integer, intent(out) :: ierr
    real(kind=REALTYPE),  intent(in)  :: y(nxl,nyl)
    real(kind=REALTYPE),  intent(out) :: yrms
    real(kind=REALTYPE) :: lsum, gsum

    ! internals
    lsum = sum(y**2)
    call MPI_Allreduce(lsum, gsum, 1, MPI_DOUBLE_PRECISION, &
                       MPI_SUM, comm, ierr)
    if (ierr /= MPI_SUCCESS) then
       write(0,*) "Error in MPI_Allreduce = ", ierr
       call MPI_Finalize(ierr)
    endif
    yrms = sqrt(gsum/nx/ny)

    ierr = 0     ! return with success flag
    return
  end subroutine PRMS
  !---------------------------------------------------------------

end module UserData
!-----------------------------------------------------------------


!-----------------------------------------------------------------
! Main driver program
!-----------------------------------------------------------------
program driver

  ! inclusions
  use UserData
  implicit none
  include "mpif.h"

  ! Declarations
  ! general problem parameters
  real*8,    parameter :: pi = 3.1415926535897932d0
  integer,   parameter :: Nt = 20           ! total number of output times 
  integer*8, parameter :: nx_ = 60          ! spatial mesh size
  integer*8, parameter :: ny_ = 120
  integer,   parameter :: PCGpretype = 1    ! enable preconditioner
  integer,   parameter :: PCGmaxl = 20      ! max num. PCG iterations
  real(kind=REALTYPE), parameter :: T0 = 0.d0        ! initial time 
  real(kind=REALTYPE), parameter :: Tf = 0.3d0       ! final time 
  real(kind=REALTYPE), parameter :: rtol = 1.d-5     ! relative and absolute tolerances
  real(kind=REALTYPE), parameter :: atol = 1.d-10
  real(kind=REALTYPE), parameter :: kx_ = 0.5d0      ! heat conductivity coefficients
  real(kind=REALTYPE), parameter :: ky_ = 0.75d0
  real(kind=REALTYPE), parameter :: nlscoef = 1.d-7  ! nonlinear solver tolerance factor

  ! solution vector and other local variables
  real(kind=REALTYPE), allocatable :: y(:,:)
  real(kind=REALTYPE) :: rout(6), rpar, t, dTout, tout, urms
  integer*8 :: iout(35), ipar, N, Ntot, i, j
  integer   :: flag, ioutput
  logical   :: outproc
  character*100 :: outname

  ! initialize MPI
  call MPI_Init(flag)
  if (flag /= MPI_SUCCESS) then
     write(0,*) "Error in MPI_Init = ", flag
     stop
  end if
  call MPI_Comm_rank(MPI_COMM_WORLD, myid, flag)
  if (flag /= MPI_SUCCESS) then
     write(0,*) "Error in MPI_Comm_rank = ", flag
     call MPI_Finalize(flag)
  end if

  ! Initialize UserData module
  call InitUserData()
  nx = nx_
  ny = ny_
  kx = kx_
  ky = ky_
  dx = 1.d0/(nx-1)   ! x mesh spacing 
  dy = 1.d0/(ny-1)   ! x mesh spacing 

  ! Set up parallel decomposition (computes local mesh sizes)
  call SetupDecomp(flag)
  if (flag /= MPI_SUCCESS) then
     write(0,*) "Error in SetupDecomp = ", flag
     call MPI_Finalize(flag)
  end if

  ! Initial problem output 
  outproc = (myid == 0)
  if (outproc) then
     write(6,*) "  "
     write(6,*) "2D Heat PDE test problem:";
     write(6,'(A,i4)') "   nprocs = " , nprocs
     write(6,'(A,i4)') "   nx = ", nx
     write(6,'(A,i4)') "   ny = ", ny
     write(6,'(A,f5.2)') "   kx = ", kx
     write(6,'(A,f5.2)') "   ky = ", ky
     write(6,'(A,es9.2)') "   rtol = ", rtol
     write(6,'(A,es9.2)') "   atol = ", atol
     write(6,'(A,i4)') "   nxl (proc 0) = ", nxl
     write(6,'(A,i4)') "   nyl (proc 0) = ", nyl
     write(6,*) "  "
  endif

  ! Initialize data structures 
  N = nxl*nyl
  Ntot = nx*ny
  call FNVInitP(comm, 4, N, Ntot, flag)
  if (flag /= MPI_SUCCESS) then
     write(0,*) "Error in FNVInitP = ", flag
     call MPI_Finalize(flag)
  end if
  allocate(y(nxl,nyl))         ! Create parallel vector for solution 
  y = 0.d0                     ! Set initial conditions 

  ! initialize PCG linear solver module
  call FSunPCGInit(4, PCGpretype, PCGmaxl, flag)
  if (flag /= MPI_SUCCESS) then
     write(0,*) "Error in FSunPCGInit = ", flag
     call MPI_Finalize(flag)
  end if
  
  ! Create the solver memory to use DIRK integrator, scalar tolerances
  call FARKMalloc(T0, y, 0, 1, rtol, atol, iout, rout, ipar, rpar, flag)
  if (flag /= MPI_SUCCESS) then
     write(0,*) "Error in FARKMalloc = ", flag
     call MPI_Finalize(flag)
  end if

  ! fill in the heat source array
  do j=1,nyl
     do i=1,nxl
        h(i,j) = sin(pi*(is+i-2)*dx) * sin(2.d0*pi*(js+j-2)*dy)
     enddo
  enddo

  ! set integrator options
  call FARKSetRin("NLCONV_COEF", nlscoef, flag)
  if (flag < 0) then
     write(0,*) "Error in FARKSetRin = ", flag
     call MPI_Finalize(flag)
  end if
  call FARKSetIin("PREDICT_METHOD", 1_8, flag)
  if (flag < 0) then
     write(0,*) "Error in FARKSetIin = ", flag
     call MPI_Finalize(flag)
  end if


  ! attach linear solver module to ARKLs interface
  call FARKLsInit(flag)
  if (flag < 0) then
     write(0,*) "Error in FARKLsInit = ", flag
     call MPI_Finalize(flag)
  end if
  call FARKLsSetPrec(1, flag)     ! Signal user-supplied preconditioner
  if (flag < 0) then
     write(0,*) "Error in FARKLsSetPrec = ", flag
     call MPI_Finalize(flag)
  end if

  ! specify that the problem is linearly implicit, but that Jacobian does not depend on time
  call FARKSetIin("LINEAR", 0_8, flag)
  if (flag < 0) then
     write(0,*) "Error in FARKSetIin = ", flag
     call MPI_Finalize(flag)
  end if

  ! Each processor outputs subdomain information
  write(outname,'(16Hheat2d_subdomain,f4.3,4H.txt)') myid/1000.0
  open(100, file=outname)
  write(100,'(6(i9,1x))') nx, ny, is, ie, js, je
  close(100)

  ! Open output streams for results, access data array 
  write(outname,'(6Hheat2d,f4.3,4H.txt)') myid/1000.0
  open(101, file=outname)

  ! Output initial condition to disk
  do j=1,nyl
     do i=1,nxl
        write(101,'(es25.16)',advance='no') y(i,j)
     enddo
  enddo
  write(101,*) "  "

  ! Main time-stepping loop: calls ARKode to perform the integration, then
  ! prints results.  Stops when the final time has been reached
  t = T0
  dTout = (Tf-T0)/Nt
  tout = T0+dTout
  call PRMS(y, urms, flag)
  if (outproc) then
    write(6,*) "        t      ||u||_rms"
    write(6,*) "   ----------------------"
    write(6,'(2(2x,f10.6))') t, urms
  endif
  do ioutput=1,Nt

     call FARKode(tout, t, y, 1, flag)         ! call integrator 
     if (flag < 0) then
        write(0,*) "Error in FARKode = ", flag
        exit
     end if
     
     call PRMS(y, urms, flag)
     if (outproc) &
          write(6,'(2(2x,f10.6))') t, urms     ! print solution stats 
     if (flag >= 0) then                       ! successful solve: update output time
        tout = min(tout + dTout, Tf)
     else                                      ! unsuccessful solve: break 
        if (outproc) &
             write(0,*) "Solver failure, stopping integration"
        exit
     endif

     ! output results to disk 
     do j=1,nyl
        do i=1,nxl
           write(101,'(es25.16)',advance='no') y(i,j)
        enddo
     enddo
     write(101,*) "  "
     
  enddo
  if (outproc) then
     write(6,*) "   ----------------------"
  endif
  close(101)

  ! Print some final statistics 
  if (outproc) then
     write(6,*) "  "
     write(6,*) "Final Solver Statistics:"
     write(6,'(2(A,i6),A)') "   Internal solver steps = ", iout(3), &
          " (attempted = ", iout(6), ")"
     write(6,'(A,i6,A,i6)') "   Total RHS evals:  Fe = ", iout(7), ",  Fi = ", iout(8)
     write(6,'(A,i6)') "   Total linear solver setups = ", iout(9)
     write(6,'(A,i6)') "   Total linear iterations = ", iout(23)
     write(6,'(A,i6)') "   Total number of Jacobian-vector products = ", iout(20)
     write(6,'(A,i6)') "   Total number of Preconditioner setups = ", iout(21)
     write(6,'(A,i6)') "   Total number of Preconditioner solves = ", iout(22)
     write(6,'(A,i6)') "   Total number of linear solver convergence failures = ", &
          iout(24)
     write(6,'(A,i6)') "   Total number of Newton iterations = ", iout(11)
     write(6,'(A,i6)') "   Total number of nonlinear solver convergence failures = ", &
          iout(12)
     write(6,'(A,i6)') "   Total number of error test failures = ", iout(10)
 endif

 ! Clean up and return with successful completion 
 if (allocated(y))  deallocate(y)    ! free solution
 call FreeUserData(flag)             ! free user data 
 call FARKFree()                     ! free integrator memory 
 call MPI_Barrier(comm, flag)
 call MPI_Finalize(flag)             ! Finalize MPI

end program driver
!-----------------------------------------------------------------



!-----------------------------------------------------------------
! Functions called by the solver
!-----------------------------------------------------------------


subroutine farkifun(t, y, ydot, ipar, rpar, ierr)
!-----------------------------------------------------------------
! f routine to compute the ODE RHS function f(t,y). 
!-----------------------------------------------------------------
  ! declarations
  use UserData
  implicit none
  include "mpif.h"
  real(kind=REALTYPE), intent(in)  :: t, rpar
  real(kind=REALTYPE), intent(in)  :: y(nxl,nyl)
  real(kind=REALTYPE), intent(out) :: ydot(nxl,nyl)
  integer*8, intent(in) :: ipar
  real(kind=REALTYPE) :: c1, c2, c3
  integer, intent(out) :: ierr
  integer*8 :: i, j
  
  ! internals

  ! Initialize ydot to zero 
  ydot = 0.d0

  ! Exchange boundary data with neighbors
  call Exchange(y, ierr)
  if (ierr /= MPI_SUCCESS) then
     write(0,*) "Error in Exchange = " , ierr
     return
  end if

  ! iterate over subdomain interior, computing approximation to RHS
  c1 = kx/dx/dx
  c2 = ky/dy/dy
  c3 = -2.d0*(c1 + c2)
  do j=2,nyl-1
     do i=2,nxl-1
        ydot(i,j) = c1*(y(i-1,j)+y(i+1,j)) + c2*(y(i,j-1)+y(i,j+1)) + c3*y(i,j)
     enddo
  enddo
  
  ! iterate over subdomain boundaries (if not at overall domain boundary)
  if (HaveNbor(1,1)) then    ! West face
     i=1
     do j=2,nyl-1
        ydot(i,j) = c1*(Wrecv(j)+y(i+1,j)) + c2*(y(i,j-1)+y(i,j+1)) + c3*y(i,j)
     enddo
  endif
  if (HaveNbor(1,2)) then    ! East face
     i=nxl
     do j=2,nyl-1
        ydot(i,j) = c1*(y(i-1,j)+Erecv(j)) + c2*(y(i,j-1)+y(i,j+1)) + c3*y(i,j)
     enddo
  endif
  if (HaveNbor(2,1)) then    ! South face
     j=1
     do i=2,nxl-1
        ydot(i,j) = c1*(y(i-1,j)+y(i+1,j)) + c2*(Srecv(i)+y(i,j+1)) + c3*y(i,j)
     enddo
  endif
  if (HaveNbor(2,2)) then    ! West face
     j=nyl
     do i=2,nxl-1
        ydot(i,j) = c1*(y(i-1,j)+y(i+1,j)) + c2*(y(i,j-1)+Nrecv(i)) + c3*y(i,j)
     enddo
  endif
  if (HaveNbor(1,1) .and. HaveNbor(2,1)) then  ! South-West corner
     i=1
     j=1
     ydot(i,j) = c1*(Wrecv(j)+y(i+1,j)) + c2*(Srecv(i)+y(i,j+1)) + c3*y(i,j)
  endif
  if (HaveNbor(1,1) .and. HaveNbor(2,2)) then  ! North-West corner
     i=1
     j=nyl
     ydot(i,j) = c1*(Wrecv(j)+y(i+1,j)) + c2*(y(i,j-1)+Nrecv(i)) + c3*y(i,j)
  endif
  if (HaveNbor(1,2) .and. HaveNbor(2,1)) then  ! South-East corner
     i=nxl 
     j=1
     ydot(i,j) = c1*(y(i-1,j)+Erecv(j)) + c2*(Srecv(i)+y(i,j+1)) + c3*y(i,j)
  endif
  if (HaveNbor(1,2) .and. HaveNbor(2,2)) then  ! North-East corner
     i=nxl
     j=nyl
     ydot(i,j) = c1*(y(i-1,j)+Erecv(j)) + c2*(y(i,j-1)+Nrecv(i)) + c3*y(i,j)
  endif

  ydot = ydot + h         ! add in heat source

  ierr = 0                ! Return with success 
  return
end subroutine farkifun
!-----------------------------------------------------------------


subroutine farkefun(t, y, ydot, ipar, rpar, ierr)
!-----------------------------------------------------------------
! Explicit portion of f routine (empty)
!-----------------------------------------------------------------
  ! declarations
  use UserData
  implicit none
  real(kind=REALTYPE), intent(in)  :: t, rpar
  integer*8, intent(in) :: ipar
  real(kind=REALTYPE), intent(in)  :: y(nxl,nyl)
  real(kind=REALTYPE), intent(out) :: ydot(nxl,nyl)
  integer, intent(out) :: ierr
  
  ! internals

  ! Initialize ydot to zero and return with success
  ydot = 0.d0
  ierr = 0
  return
end subroutine farkefun
!-----------------------------------------------------------------


subroutine farkpset(t, y, fy, jok, jcur, gamma, hcur, ipar, &
                    rpar, ierr)
!-----------------------------------------------------------------
! Preconditioner setup routine (fills inverse of Jacobian diagonal)
!-----------------------------------------------------------------
  ! declarations
  use UserData
  implicit none
  real(kind=REALTYPE), intent(in) :: t, gamma, hcur, rpar
  real(kind=REALTYPE), intent(in) :: y(nxl,nyl), fy(nxl,nyl)
  integer*8, intent(in) :: ipar
  integer, intent(in) :: jok
  integer, intent(out) :: jcur, ierr
  real(kind=REALTYPE) :: c

  ! internals
  c = 1.d0 + gamma*2.d0*(kx/dx/dx + ky/dy/dy)

  ! set all entries of d to the inverse of the diagonal values in interior
  ! (since boundary RHS is 0, set boundary diagonals to the same)
  d = 1.d0/c

  jcur = 1     ! update jcur flag
  ierr = 0     ! Return with success 
  return
end subroutine farkpset


subroutine farkpsol(t, y, fy, r, z, gamma, delta, lr, &
                    ipar, rpar, ierr)
!-----------------------------------------------------------------
! Preconditioner solve routine
!-----------------------------------------------------------------
  ! declarations
  use UserData
  implicit none
  real(kind=REALTYPE), intent(in)  :: t, gamma, delta, rpar
  integer*8, intent(in) :: ipar
  real(kind=REALTYPE), intent(in)  :: y(nxl,nyl), fy(nxl,nyl), r(nxl,nyl)
  real(kind=REALTYPE), intent(out) :: z(nxl,nyl)
  integer, intent(in)  :: lr
  integer, intent(out) :: ierr

  ! internals
  z = r*d      ! perform Jacobi iteration (whole array operation)
  ierr = 0     ! Return with success 
  return
end subroutine farkpsol
!-----------------------------------------------------------------

