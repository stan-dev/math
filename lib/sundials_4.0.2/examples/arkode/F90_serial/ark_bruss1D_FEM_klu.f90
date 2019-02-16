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
! The following test simulates a brusselator problem from chemical 
! kinetics.  This is a PDE system with 3 components, Y = [u,v,w], 
! satisfying the equations,
!    du/dt = du*u_xx + a - (w+1)*u + v*u^2
!    dv/dt = dv*v_xx + w*u - v*u^2
!    dw/dt = dw*w_xx + (b-w)/ep - w*u
! for t in the interval [0, 10], x in [0, 10], with initial 
! conditions 
!    u(0,x) =  a  + 0.1*sin(pi*x),
!    v(0,x) = b/a + 0.1*sin(pi*x),
!    w(0,x) =  b  + 0.1*sin(pi*x),
! and with stationary boundary conditions, i.e.
!    u_t(t,0) = u_t(t,1) = 0,
!    v_t(t,0) = v_t(t,1) = 0,
!    w_t(t,0) = w_t(t,1) = 0.
!
! Here, we use a piecewise linear Galerkin finite element 
! discretization in space, where all element-wise integrals are 
! computed using 3-node Gaussian quadrature (since we will have 
! quartic polynomials in the reaction terms for the u_t and v_t 
! equations (including the test function)).  The time derivative 
! terms in this system will include a mass matrix, giving rise to 
! an ODE system of the form
!      M y_t = L y + R(y),
! where M is the 3x3 block mass matrix for each component, L is 
! the 3x3 block Laplace operator for each component, and R(y) is 
! comprised of the nonlinear reaction terms for each component.  
! Since it it highly inefficient to rewrite this system as
!      y_t = M^{-1}(L y + R(y)),
! we solve this system using FARKODE, with a user-supplied mass
! matrix.  We therefore provide functions to evaluate the ODE RHS 
!    f(t,y) = L y + R(y),
! its Jacobian
!    J(t,y) = L + dR/dy,
! and the mass matrix, M.
!
! We use N=201 spatial nodes, with parameters 
!    a=0.6,  b=2.0,  du=0.025,  dv=0.025,  dw=0.025,  ep=1.d-5
!
! This program solves the problem with the DIRK method, using a 
! Newton iteration with the SUNKLU sparse linear solvers for both 
! the system and mass matrices.  These matrices are stored in 
! compressed-sparse-row format.
!
! Output is printed 10 times throughout the defined time interval.
! Run statistics (optional outputs) are printed at the end.
!-----------------------------------------------------------------


! user data structure
module UserData
  implicit none
  include "sundials/sundials_fconfig.h"
  save
  
  integer*8 :: N                   ! number of intervals
  real*8, allocatable :: x(:)      ! mesh node locations
  real*8 :: a                      ! constant forcing on u
  real*8 :: b                      ! steady-state value of w
  real*8 :: du                     ! diffusion coeff for u
  real*8 :: dv                     ! diffusion coeff for v
  real*8 :: dw                     ! diffusion coeff for w
  real*8 :: ep                     ! stiffness parameter

contains

  ! function that maps 2D data into 1D address space
  ! (0-based since CSR matrix will be sent to C solver)
  integer(kind=8) function idx(ix,ivar)
    integer :: ivar, ix
    idx = 3*(ix-1) + ivar - 1
  end function idx

end module UserData


! finite element basis functions
module FEM
  
contains
  
  ! left/right basis functions
  double precision function ChiL(xl,xr,x)
    double precision :: xl, xr, x
    ChiL = (xr-x)/(xr-xl)
  end function ChiL

  double precision function ChiR(xl,xr,x)
    double precision :: xl, xr, x
    ChiR = (x-xl)/(xr-xl)
  end function ChiR

  ! derivatives of left/right basis functions
  double precision function ChiL_x(xl,xr)
    double precision :: xl, xr
    ChiL_x = 1.d0/(xl-xr)
  end function ChiL_X

  double precision function ChiR_x(xl,xr)
    double precision :: xl, xr
    ChiR_x = 1.d0/(xr-xl)
  end function ChiR_x

  ! FEM output evaluation routines: value and derivative
  double precision function Eval(ul,ur,xl,xr,x)
    double precision :: ul, ur, xl, xr, x
    Eval = ul*ChiL(xl,xr,x) + ur*ChiR(xl,xr,x)
  end function Eval

  double precision function Eval_x(ul,ur,xl,xr)
    double precision :: ul, ur, xl, xr
    Eval_x = ul*ChiL_x(xl,xr) + ur*ChiR_x(xl,xr)
  end function Eval_x

end module FEM


! quadrature data
module Quadrature
  
contains

  ! nodes
  double precision function X1(xl,xr)
    double precision :: xl, xr
    X1 = 0.5d0*(xl+xr) - 0.5d0*(xr-xl)*0.774596669241483377035853079956d0
  end function X1

  double precision function X2(xl,xr)
    double precision :: xl, xr
    X2 = 0.5d0*(xl+xr)
  end function X2

  double precision function X3(xl,xr)
    double precision :: xl, xr
    X3 = 0.5d0*(xl+xr) + 0.5d0*(xr-xl)*0.774596669241483377035853079956d0
  end function X3

  ! quadrature 
  double precision function Quad(f1,f2,f3,xl,xr)
    real*8 :: f1, f2, f3, xl, xr
    real*8, parameter :: wt1=0.55555555555555555555555555555556d0
    real*8, parameter :: wt2=0.88888888888888888888888888888889d0
    real*8, parameter :: wt3=0.55555555555555555555555555555556d0
    Quad = 0.5d0*(xr-xl)*(wt1*f1 + wt2*f2 + wt3*f3)
  end function Quad

end module Quadrature
!-----------------------------------------------------------------




!-----------------------------------------------------------------
! Main driver program
!-----------------------------------------------------------------
program driver

  ! inclusions
  use UserData
  implicit none

  ! Declarations
  ! general problem variables
  real(kind=REALTYPE), parameter :: T0=0.d0, Tf=10.d0
  real(kind=REALTYPE) :: rtol, atol, rout(6), Tout, Tcur
  real*8    :: dTout, pi, h, z
  integer   :: i, it, Nt, ier, ordering, sparsetype, time_dep
  integer*8 :: iout(35), NEQ, nnz, Iinput
  real(kind=REALTYPE), allocatable :: y(:,:), umask(:,:)
  real(kind=REALTYPE), allocatable :: vmask(:,:), wmask(:,:)

  ! dummy real/integer parameters to pass through to supplied functions
  integer*8 :: ipar
  real(kind=REALTYPE) :: rpar

  !-----------------------

  ! set problem parameters (in UserData structure)
  N = 201                 ! number of intervals
  a = 0.6d0               ! constant forcing on u
  b = 2.d0                ! steady-state value of w
  du = 2.5d-2             ! diffusion coeff for u
  dv = 2.5d-2             ! diffusion coeff for v
  dw = 2.5d-2             ! diffusion coeff for w
  ep = 1.d-5              ! stiffness parameter

  ! set overall problem size, allocate solution/mask/temporary arrays
  NEQ = 3*N
  allocate(y(3,N), umask(3,N), vmask(3,N), wmask(3,N))

  ! allocate and set up spatial mesh; this [arbitrarily] clusters 
  ! more intervals near the end points of the interval
  allocate(x(N))          ! mesh node locations
  pi = 4.d0*atan(1.d0)
  h = 10.d0/(N-1)
  do i=1,N
     z = -5.d0 + h*(i-1)
     x(i) = 0.5d0/atan(5.d0)*atan(z) + 0.5d0
  enddo

  ! output mesh to disk
  open(200, file='bruss_FEM_mesh.txt')
  do i=1,N
     write(200,*) x(i)
  enddo
  close(200)

  ! time-stepping information
  dTout = (Tf-T0)/10.d0
  Nt = Tf/dTout + 0.5

  ! set initial conditions
  do i=1,N
     y(1,i) =  a  + 0.1d0*sin(pi*x(i))   ! u0
     y(2,i) = b/a + 0.1d0*sin(pi*x(i))   ! v0
     y(3,i) =  b  + 0.1d0*sin(pi*x(i))   ! w0
  enddo

  ! set mask values for each solution component
  umask = 0.d0
  vmask = 0.d0
  wmask = 0.d0
  do i=1,N
     umask(1,i) = 1.d0
     vmask(2,i) = 1.d0
     wmask(3,i) = 1.d0
  enddo

  ! set tolerances according to problem specifications
  atol = 1.d-11
  rtol = 1.d-6
  
  ! initialize vector module
  call FNVInitS(4, NEQ, ier)

  ! initialize system and mass matrix modules
  nnz = 15*NEQ     ! integer number of nonzeros           
  ordering = 0     ! AMD
  sparsetype = 1   ! CSR
  call FSunSparseMatInit(4, NEQ, NEQ, nnz, sparsetype, ier)
  call FSunSparseMassMatInit(NEQ, NEQ, nnz, sparsetype, ier)

  ! initialize KLU system and mass solvers
  call FSunKLUInit(4, ier)
  call FSunMassKLUInit(ier)
  
  ! initialize ARKode solver
  ipar = 0
  rpar = 0.0
  iout = 0
  rout = 0.0
  call FARKMalloc(T0, y, 0, 1, rtol, atol, &
                  iout, rout, ipar, rpar, ier)

  ! set optional inputs
  Iinput = 1
  call FARKSetIin('IMPLICIT', Iinput, ier)
  Iinput = 1000
  call FARKSetIin('MAX_NSTEPS', Iinput, ier)
  call FARKSetResTolerance(1, atol, ier)

  ! attach matrix and linear solver objects to ARKLs interfaces
  time_dep = 0
  call FARKLsInit(ier)
  call FARKSparseSetJac(ier)
  call FARKLsMassInit(time_dep, ier)
  call FARKSparseSetMass(ier)

  ! Open output stream for results
  open(501, file='bruss_FEM_u.txt')
  open(502, file='bruss_FEM_v.txt')
  open(503, file='bruss_FEM_w.txt')

  ! output initial condition to disk
  write(501,*) ( y(1,i), i=1,N )
  write(502,*) ( y(2,i), i=1,N )
  write(503,*) ( y(3,i), i=1,N )
 
  ! output solver parameters to screen
  call FARKWriteParameters(ier)

  ! loop over time outputs
  Tout = T0
  Tcur = T0
  print *, '        t         ||u||_rms    ||v||_rms    ||w||_rms'
  print *, '  ----------------------------------------------------'
  print '(3x,4(es12.5,1x))', Tcur, sqrt(sum(y*y*umask)/N), &
       sqrt(sum(y*y*vmask)/N), sqrt(sum(y*y*wmask)/N)
  do it = 1,Nt

     ! set next output time
     Tout = min(Tout + dTout, Tf)

     ! set next output time
     call FARKSetRin('STOP_TIME', Tout, ier)

     ! call solver
     call FARKode(Tout, Tcur, y, 1, ier)
     if (ier < 0) then
        write(0,*) 'Solver failure, stopping integration'
        stop
     end if

     ! output current solution information
     print '(3x,4(es12.5,1x))', Tcur, sqrt(sum(y*y*umask)/N), &
          sqrt(sum(y*y*vmask)/N), sqrt(sum(y*y*wmask)/N)

     ! output current results to disk
     write(501,*) ( y(1,i), i=1,N )
     write(502,*) ( y(2,i), i=1,N )
     write(503,*) ( y(3,i), i=1,N )
 
  end do
  print *, '  ----------------------------------------------------'

  ! close solution output files
  close(501)
  close(502)
  close(503)

  ! output solver statistics
  print *, '  '
  print *, 'Final Solver Statistics:'
  print *, '   Internal solver steps = ', iout(3),' (attempted =', iout(6), ')'
  print *, '   Total RHS evals:  Fe = ', iout(7),'  Fi = ', iout(8)
  print *, '   Total linear solver setups = ', iout(9)
  print *, '   Total RHS evals for setting up the linear system = ', iout(17)
  print *, '   Total number of Jacobian evaluations = ', iout(18)
  print *, '   Total number of nonlinear iterations = ', iout(11)
  print *, '   Total number of nonlinear solver convergence failures = ',iout(12)
  print *, '   Total number of error test failures = ', iout(10)
  print *, '   Total number of mass matrix evaluations = ', iout(28)
  print *, '  '

  ! clean up
  deallocate(y, umask, vmask, wmask, x)
  call FARKFree()

end program driver
!-----------------------------------------------------------------




!-----------------------------------------------------------------
! User-supplied subroutines for FARKODE interface
!-----------------------------------------------------------------


subroutine FARKIFun(t, y, ydot, ipar, rpar, ier)
!-----------------------------------------------------------------
! Right-hand side of the ODE system
!-----------------------------------------------------------------
  use UserData
  use FEM
  use Quadrature

  ! Declarations
  implicit none

  ! Arguments
  real(kind=REALTYPE), intent(in)  :: t, rpar(1)
  integer*8, intent(in) :: ipar(1)
  integer, intent(out)  :: ier
  real(kind=REALTYPE), intent(in)  :: y(3,N)
  real(kind=REALTYPE), intent(out) :: ydot(3,N)

  ! Local data
  integer :: ix
  logical :: left, right
  real*8  :: ul, ur, vl, vr, wl, wr, xl, xr, u, v, w, f1, f2, f3

  ! clear out rhs
  ydot = 0.d0

  ! iterate over intervals, filling in rhs function
  do ix=1,N-1

     ! set booleans to determine whether equations exist on the left/right */
     left  = .true.
     right = .true.
     if (ix==1)      left  = .false.
     if (ix==(N-1))  right = .false.

     ! set nodal value shortcuts (interval index aligns with left node)
     ul = y(1,ix)
     vl = y(2,ix)
     wl = y(3,ix)
     ur = y(1,ix+1)
     vr = y(2,ix+1)
     wr = y(3,ix+1)

     ! set mesh shortcuts
     xl = x(ix)
     xr = x(ix+1)

     !    left test function
     if (left) then

        ! u -- reaction
        u = Eval(ul, ur, xl, xr, X1(xl,xr))
        v = Eval(vl, vr, xl, xr, X1(xl,xr))
        w = Eval(wl, wr, xl, xr, X1(xl,xr))
        f1 = (a - (w+1.d0)*u + v*u*u) * ChiL(xl,xr,X1(xl,xr))
        u = Eval(ul, ur, xl, xr, X2(xl,xr))
        v = Eval(vl, vr, xl, xr, X2(xl,xr))
        w = Eval(wl, wr, xl, xr, X2(xl,xr))
        f2 = (a - (w+1.d0)*u + v*u*u) * ChiL(xl,xr,X2(xl,xr))
        u = Eval(ul, ur, xl, xr, X3(xl,xr))
        v = Eval(vl, vr, xl, xr, X3(xl,xr))
        w = Eval(wl, wr, xl, xr, X3(xl,xr))
        f3 = (a - (w+1.d0)*u + v*u*u) * ChiL(xl,xr,X3(xl,xr))
        ydot(1,ix) = ydot(1,ix) + Quad(f1,f2,f3,xl,xr)

        ! u -- diffusion
        f1 = -du * Eval_x(ul,ur,xl,xr) * ChiL_x(xl,xr)
        ydot(1,ix) = ydot(1,ix) + Quad(f1,f1,f1,xl,xr)

        ! v -- reaction
        u = Eval(ul, ur, xl, xr, X1(xl,xr))
        v = Eval(vl, vr, xl, xr, X1(xl,xr))
        w = Eval(wl, wr, xl, xr, X1(xl,xr))
        f1 = (w*u - v*u*u) * ChiL(xl,xr,X1(xl,xr))
        u = Eval(ul, ur, xl, xr, X2(xl,xr))
        v = Eval(vl, vr, xl, xr, X2(xl,xr))
        w = Eval(wl, wr, xl, xr, X2(xl,xr))
        f2 = (w*u - v*u*u) * ChiL(xl,xr,X2(xl,xr))
        u = Eval(ul, ur, xl, xr, X3(xl,xr))
        v = Eval(vl, vr, xl, xr, X3(xl,xr))
        w = Eval(wl, wr, xl, xr, X3(xl,xr))
        f3 = (w*u - v*u*u) * ChiL(xl,xr,X3(xl,xr))
        ydot(2,ix) = ydot(2,ix) + Quad(f1,f2,f3,xl,xr)
      
        ! v -- diffusion
        f1 = -dv * Eval_x(vl,vr,xl,xr) * ChiL_x(xl,xr)
        ydot(2,ix) = ydot(2,ix) + Quad(f1,f1,f1,xl,xr)
      
        ! w -- reaction
        u = Eval(ul, ur, xl, xr, X1(xl,xr))
        v = Eval(vl, vr, xl, xr, X1(xl,xr))
        w = Eval(wl, wr, xl, xr, X1(xl,xr))
        f1 = ((b-w)/ep - w*u) * ChiL(xl,xr,X1(xl,xr))
        u = Eval(ul, ur, xl, xr, X2(xl,xr))
        v = Eval(vl, vr, xl, xr, X2(xl,xr))
        w = Eval(wl, wr, xl, xr, X2(xl,xr))
        f2 = ((b-w)/ep - w*u) * ChiL(xl,xr,X2(xl,xr))
        u = Eval(ul, ur, xl, xr, X3(xl,xr))
        v = Eval(vl, vr, xl, xr, X3(xl,xr))
        w = Eval(wl, wr, xl, xr, X3(xl,xr))
        f3 = ((b-w)/ep - w*u) * ChiL(xl,xr,X3(xl,xr))
        ydot(3,ix) = ydot(3,ix) + Quad(f1,f2,f3,xl,xr)

        ! w -- diffusion
        f1 = -dw * Eval_x(wl,wr,xl,xr) * ChiL_x(xl,xr)
        ydot(3,ix) = ydot(3,ix) + Quad(f1,f1,f1,xl,xr)

     end if

     !    right test function 
     if (right) then

        ! u -- reaction
        u = Eval(ul, ur, xl, xr, X1(xl,xr))
        v = Eval(vl, vr, xl, xr, X1(xl,xr))
        w = Eval(wl, wr, xl, xr, X1(xl,xr))
        f1 = (a - (w+1.d0)*u + v*u*u) * ChiR(xl,xr,X1(xl,xr))
        u = Eval(ul, ur, xl, xr, X2(xl,xr))
        v = Eval(vl, vr, xl, xr, X2(xl,xr))
        w = Eval(wl, wr, xl, xr, X2(xl,xr))
        f2 = (a - (w+1.d0)*u + v*u*u) * ChiR(xl,xr,X2(xl,xr))
        u = Eval(ul, ur, xl, xr, X3(xl,xr))
        v = Eval(vl, vr, xl, xr, X3(xl,xr))
        w = Eval(wl, wr, xl, xr, X3(xl,xr))
        f3 = (a - (w+1.d0)*u + v*u*u) * ChiR(xl,xr,X3(xl,xr))
        ydot(1,ix+1) = ydot(1,ix+1) + Quad(f1,f2,f3,xl,xr)

        ! u -- diffusion
        f1 = -du * Eval_x(ul,ur,xl,xr) * ChiR_x(xl,xr)
        ydot(1,ix+1) = ydot(1,ix+1) + Quad(f1,f1,f1,xl,xr)

        ! v -- reaction
        u = Eval(ul, ur, xl, xr, X1(xl,xr))
        v = Eval(vl, vr, xl, xr, X1(xl,xr))
        w = Eval(wl, wr, xl, xr, X1(xl,xr))
        f1 = (w*u - v*u*u) * ChiR(xl,xr,X1(xl,xr))
        u = Eval(ul, ur, xl, xr, X2(xl,xr))
        v = Eval(vl, vr, xl, xr, X2(xl,xr))
        w = Eval(wl, wr, xl, xr, X2(xl,xr))
        f2 = (w*u - v*u*u) * ChiR(xl,xr,X2(xl,xr))
        u = Eval(ul, ur, xl, xr, X3(xl,xr))
        v = Eval(vl, vr, xl, xr, X3(xl,xr))
        w = Eval(wl, wr, xl, xr, X3(xl,xr))
        f3 = (w*u - v*u*u) * ChiR(xl,xr,X3(xl,xr))
        ydot(2,ix+1) = ydot(2,ix+1) + Quad(f1,f2,f3,xl,xr)

        ! v -- diffusion
        f1 = -dv * Eval_x(vl,vr,xl,xr) * ChiR_x(xl,xr)
        ydot(2,ix+1) = ydot(2,ix+1) + Quad(f1,f1,f1,xl,xr)

        ! w -- reaction
        u = Eval(ul, ur, xl, xr, X1(xl,xr))
        v = Eval(vl, vr, xl, xr, X1(xl,xr))
        w = Eval(wl, wr, xl, xr, X1(xl,xr))
        f1 = ((b-w)/ep - w*u) * ChiR(xl,xr,X1(xl,xr))
        u = Eval(ul, ur, xl, xr, X2(xl,xr))
        v = Eval(vl, vr, xl, xr, X2(xl,xr))
        w = Eval(wl, wr, xl, xr, X2(xl,xr))
        f2 = ((b-w)/ep - w*u) * ChiR(xl,xr,X2(xl,xr))
        u = Eval(ul, ur, xl, xr, X3(xl,xr))
        v = Eval(vl, vr, xl, xr, X3(xl,xr))
        w = Eval(wl, wr, xl, xr, X3(xl,xr))
        f3 = ((b-w)/ep - w*u) * ChiR(xl,xr,X3(xl,xr))
        ydot(3,ix+1) = ydot(3,ix+1) + Quad(f1,f2,f3,xl,xr)

        ! w -- diffusion
        f1 = -dw * Eval_x(wl,wr,xl,xr) * ChiR_x(xl,xr)
        ydot(3,ix+1) = ydot(3,ix+1) + Quad(f1,f1,f1,xl,xr)

     endif

  enddo
  ier = 0

  return

end subroutine farkifun
!-----------------------------------------------------------------



subroutine farkefun(t, y, ydot, ipar, rpar, ier)
!-----------------------------------------------------------------
! (unused) Explicit portion of the ODE right-hand function
!-----------------------------------------------------------------
  use UserData

  ! Declarations
  implicit none

  ! Arguments
  real(kind=REALTYPE), intent(in)  :: t, rpar(1)
  integer*8, intent(in) :: ipar(1)
  real(kind=REALTYPE), intent(in)  :: y(3,N)
  real(kind=REALTYPE), intent(out) :: ydot(3,N)
  integer, intent(out)  :: ier

  ! return with success (since fully implicit)
  ydot = 0.d0
  ier = 0
  
end subroutine farkefun
!-----------------------------------------------------------------



subroutine farkspjac(t, y, fy, neq, nnz, Jdata, Jcolvals, &
                     Jrowptrs, h, ipar, rpar, wk1, wk2, wk3, ier)
!-----------------------------------------------------------------
! Jacobian computation routine -- CSR format
!-----------------------------------------------------------------
  use UserData
  use Quadrature
  use FEM

  ! Declarations
  implicit none

  ! Arguments
  real(kind=REALTYPE), intent(in)  :: t, h, rpar(1)
  real(kind=REALTYPE), intent(in), dimension(3,N) :: y, fy, wk1, wk2, wk3
  real(kind=REALTYPE), intent(out) :: Jdata(nnz)
  integer*8, intent(in)  :: ipar(1), neq, nnz
  integer(kind=SUNINDEXTYPE), intent(out) :: Jcolvals(nnz)
  integer(kind=SUNINDEXTYPE), intent(out) :: Jrowptrs(neq+1)
  integer,   intent(out) :: ier

  ! Local data
  integer :: ix, nz, Nint
  real*8  :: ul, uc, ur, vl, vc, vr, wl, wc, wr, xl, xc, xr
  real*8  :: u1, u2, u3, v1, v2, v3, w1, w2, w3
  real*8  :: f1, f2, f3, df1, df2, df3, dQdf1, dQdf2, dQdf3
  real*8  :: ChiL1, ChiL2, ChiL3, ChiR1, ChiR2, ChiR3
  real*8, dimension(3,-1:1) :: Ju, Jv, Jw

  ! check that vector/matrix dimensions match up
  if ((3*N /= neq) .or. (nnz < 27*(N-2))) then
     ier = 1
     return
  endif

  ! set integer*4 version of N for call to idx()
  Nint = N
  
  ! clear out Jacobian matrix data
  Jdata = 0.d0
  nz = 0

  ! Dirichlet boundary at left
  Jrowptrs(idx(1,1)+1) = nz
  Jrowptrs(idx(1,2)+1) = nz
  Jrowptrs(idx(1,3)+1) = nz
 
  ! iterate through nodes, filling in matrix by rows
  do ix=2,N-1

     ! set nodal value shortcuts (interval index aligns with left node)
     xl = x(ix-1)
     ul = y(1,ix-1)
     vl = y(2,ix-1)
     wl = y(3,ix-1)
     xc = x(ix)
     uc = y(1,ix)
     vc = y(2,ix)
     wc = y(3,ix)
     xr = x(ix+1)
     ur = y(1,ix+1)
     vr = y(2,ix+1)
     wr = y(3,ix+1)

     ! compute entries of all Jacobian rows at node ix
     Ju = 0.d0
     Jv = 0.d0
     Jw = 0.d0

     ! first compute dependence on values to left and center
     
     !    evaluate relevant variables in left subinterval
     u1 = Eval(ul, uc, xl, xc, X1(xl,xc))
     v1 = Eval(vl, vc, xl, xc, X1(xl,xc))
     w1 = Eval(wl, wc, xl, xc, X1(xl,xc))
     u2 = Eval(ul, uc, xl, xc, X2(xl,xc))
     v2 = Eval(vl, vc, xl, xc, X2(xl,xc))
     w2 = Eval(wl, wc, xl, xc, X2(xl,xc))
     u3 = Eval(ul, uc, xl, xc, X3(xl,xc))
     v3 = Eval(vl, vc, xl, xc, X3(xl,xc))
     w3 = Eval(wl, wc, xl, xc, X3(xl,xc))

     dQdf1 = Quad(1.d0, 0.d0, 0.d0, xl, xc)
     dQdf2 = Quad(0.d0, 1.d0, 0.d0, xl, xc)
     dQdf3 = Quad(0.d0, 0.d0, 1.d0, xl, xc)

     ChiL1 = ChiL(xl, xc, X1(xl,xc))
     ChiL2 = ChiL(xl, xc, X2(xl,xc))
     ChiL3 = ChiL(xl, xc, X3(xl,xc))
     ChiR1 = ChiR(xl, xc, X1(xl,xc))
     ChiR2 = ChiR(xl, xc, X2(xl,xc))
     ChiR3 = ChiR(xl, xc, X3(xl,xc))


     !    compute diffusion Jacobian components

     !    L_u = -du * u_x * ChiR_x
     !     dL_u/dul
     Ju(1,-1) = (-du) * Quad(1.d0,1.d0,1.d0,xl,xc) * ChiL_x(xl,xc) * ChiR_x(xl,xc)
     !     dL_u/duc
     Ju(1,0)  = (-du) * Quad(1.d0,1.d0,1.d0,xl,xc) * ChiR_x(xl,xc) * ChiR_x(xl,xc)

     !    L_v = -dv * v_x * ChiR_x
     !     dL_v/dvl
     Jv(2,-1) = (-dv) * Quad(1.d0,1.d0,1.d0,xl,xc) * ChiL_x(xl,xc) * ChiR_x(xl,xc)
     !     dL_v/dvc
     Jv(2,0)  = (-dv) * Quad(1.d0,1.d0,1.d0,xl,xc) * ChiR_x(xl,xc) * ChiR_x(xl,xc)

     !    L_w =  -dw * w_x * ChiR_x
     !     dL_w/dwl
     Jw(3,-1) = (-dw) * Quad(1.d0,1.d0,1.d0,xl,xc) * ChiL_x(xl,xc) * ChiR_x(xl,xc)
     !     dL_w/dwc
     Jw(3,0)  = (-dw) * Quad(1.d0,1.d0,1.d0,xl,xc) * ChiR_x(xl,xc) * ChiR_x(xl,xc)


     !    compute reaction Jacobian components

     !    R_u = (a - (w+1.d0)*u + v*u*u) 
     !     dR_u/dul 
     df1 = (-(w1+1.d0) + 2.d0*v1*u1) * ChiL1 * ChiR1
     df2 = (-(w2+1.d0) + 2.d0*v2*u2) * ChiL2 * ChiR2
     df3 = (-(w3+1.d0) + 2.d0*v3*u3) * ChiL3 * ChiR3
     Ju(1,-1) = Ju(1,-1) + dQdf1*df1 + dQdf2*df2 + dQdf3*df3

     !     dR_u/duc
     df1 = (-(w1+1.d0) + 2.d0*v1*u1) * ChiR1 * ChiR1
     df2 = (-(w2+1.d0) + 2.d0*v2*u2) * ChiR2 * ChiR2
     df3 = (-(w3+1.d0) + 2.d0*v3*u3) * ChiR3 * ChiR3
     Ju(1,0) = Ju(1,0)+ dQdf1*df1 + dQdf2*df2 + dQdf3*df3

     !     dR_u/dvl 
     df1 = (u1*u1) * ChiL1 * ChiR1
     df2 = (u2*u2) * ChiL2 * ChiR2
     df3 = (u3*u3) * ChiL3 * ChiR3
     Ju(2,-1) = Ju(2,-1) + dQdf1*df1 + dQdf2*df2 + dQdf3*df3

     !     dR_u/dvc
     df1 = (u1*u1) * ChiR1 * ChiR1
     df2 = (u2*u2) * ChiR2 * ChiR2
     df3 = (u3*u3) * ChiR3 * ChiR3
     Ju(2,0) = Ju(2,0) + dQdf1*df1 + dQdf2*df2 + dQdf3*df3

     !     dR_u/dwl 
     df1 = (-u1) * ChiL1 * ChiR1
     df2 = (-u2) * ChiL2 * ChiR2
     df3 = (-u3) * ChiL3 * ChiR3
     Ju(3,-1) = Ju(3,-1) + dQdf1*df1 + dQdf2*df2 + dQdf3*df3

     !     dR_u/dwc
     df1 = (-u1) * ChiR1 * ChiR1
     df2 = (-u2) * ChiR2 * ChiR2
     df3 = (-u3) * ChiR3 * ChiR3
     Ju(3,0) = Ju(3,0) + dQdf1*df1 + dQdf2*df2 + dQdf3*df3


     !    R_v = (w*u - v*u*u) 
     !     dR_v/dul 
     df1 = (w1 - 2.d0*v1*u1) * ChiL1 * ChiR1
     df2 = (w2 - 2.d0*v2*u2) * ChiL2 * ChiR2
     df3 = (w3 - 2.d0*v3*u3) * ChiL3 * ChiR3
     Jv(1,-1) = Jv(1,-1) + dQdf1*df1 + dQdf2*df2 + dQdf3*df3

     !     dR_v/duc
     df1 = (w1 - 2.d0*v1*u1) * ChiR1 * ChiR1
     df2 = (w2 - 2.d0*v2*u2) * ChiR2 * ChiR2
     df3 = (w3 - 2.d0*v3*u3) * ChiR3 * ChiR3
     Jv(1,0) = Jv(1,0) + dQdf1*df1 + dQdf2*df2 + dQdf3*df3

     !     dR_v/dvl 
     df1 = (-u1*u1) * ChiL1 * ChiR1
     df2 = (-u2*u2) * ChiL2 * ChiR2
     df3 = (-u3*u3) * ChiL3 * ChiR3
     Jv(2,-1) = Jv(2,-1) + dQdf1*df1 + dQdf2*df2 + dQdf3*df3

     !     dR_v/dvc
     df1 = (-u1*u1) * ChiR1 * ChiR1
     df2 = (-u2*u2) * ChiR2 * ChiR2
     df3 = (-u3*u3) * ChiR3 * ChiR3
     Jv(2,0) = Jv(2,0) + dQdf1*df1 + dQdf2*df2 + dQdf3*df3

     !     dR_v/dwl 
     df1 = (u1) * ChiL1 * ChiR1
     df2 = (u2) * ChiL2 * ChiR2
     df3 = (u3) * ChiL3 * ChiR3
     Jv(3,-1) = Jv(3,-1) + dQdf1*df1 + dQdf2*df2 + dQdf3*df3

     !     dR_v/dwc 
     df1 = (u1) * ChiR1 * ChiR1
     df2 = (u2) * ChiR2 * ChiR2
     df3 = (u3) * ChiR3 * ChiR3
     Jv(3,0) = Jv(3,0) + dQdf1*df1 + dQdf2*df2 + dQdf3*df3


     !    R_w = ((b-w)/ep - w*u) 
     !     dR_w/dul 
     df1 = (-w1) * ChiL1 * ChiR1
     df2 = (-w2) * ChiL2 * ChiR2
     df3 = (-w3) * ChiL3 * ChiR3
     Jw(1,-1) = Jw(1,-1) + dQdf1*df1 + dQdf2*df2 + dQdf3*df3

     !     dR_w/duc
     df1 = (-w1) * ChiR1 * ChiR1
     df2 = (-w2) * ChiR2 * ChiR2
     df3 = (-w3) * ChiR3 * ChiR3
     Jw(1,0) = Jw(1,0) + dQdf1*df1 + dQdf2*df2 + dQdf3*df3

     !     dR_w/dwl 
     df1 = (-1.d0/ep - u1) * ChiL1 * ChiR1
     df2 = (-1.d0/ep - u2) * ChiL2 * ChiR2
     df3 = (-1.d0/ep - u3) * ChiL3 * ChiR3
     Jw(3,-1) = Jw(3,-1) + dQdf1*df1 + dQdf2*df2 + dQdf3*df3

     !     dR_w/dwc
     df1 = (-1.d0/ep - u1) * ChiR1 * ChiR1
     df2 = (-1.d0/ep - u2) * ChiR2 * ChiR2
     df3 = (-1.d0/ep - u3) * ChiR3 * ChiR3
     Jw(3,0) = Jw(3,0) + dQdf1*df1 + dQdf2*df2 + dQdf3*df3


     ! second compute dependence on values to center and right

     !    evaluate relevant variables in right subinterval
     u1 = Eval(uc, ur, xc, xr, X1(xc,xr))
     v1 = Eval(vc, vr, xc, xr, X1(xc,xr))
     w1 = Eval(wc, wr, xc, xr, X1(xc,xr))
     u2 = Eval(uc, ur, xc, xr, X2(xc,xr))
     v2 = Eval(vc, vr, xc, xr, X2(xc,xr))
     w2 = Eval(wc, wr, xc, xr, X2(xc,xr))
     u3 = Eval(uc, ur, xc, xr, X3(xc,xr))
     v3 = Eval(vc, vr, xc, xr, X3(xc,xr))
     w3 = Eval(wc, wr, xc, xr, X3(xc,xr))

     dQdf1 = Quad(1.d0, 0.d0, 0.d0, xc, xr)
     dQdf2 = Quad(0.d0, 1.d0, 0.d0, xc, xr)
     dQdf3 = Quad(0.d0, 0.d0, 1.d0, xc, xr)

     ChiL1 = ChiL(xc, xr, X1(xc,xr))
     ChiL2 = ChiL(xc, xr, X2(xc,xr))
     ChiL3 = ChiL(xc, xr, X3(xc,xr))
     ChiR1 = ChiR(xc, xr, X1(xc,xr))
     ChiR2 = ChiR(xc, xr, X2(xc,xr))
     ChiR3 = ChiR(xc, xr, X3(xc,xr))


     !    compute diffusion Jacobian components

     !    L_u = -du * u_x * ChiL_x
     !     dL_u/duc
     Ju(1,0) = Ju(1,0) + (-du) * Quad(1.d0,1.d0,1.d0,xc,xr) * ChiL_x(xc,xr) * ChiL_x(xc,xr)

     !     dL_u/dur
     Ju(1,1) = Ju(1,1) + (-du) * Quad(1.d0,1.d0,1.d0,xc,xr) * ChiL_x(xc,xr) * ChiR_x(xc,xr)

     !    L_v = -dv * v_x * ChiL_x
     !     dL_v/dvc
     Jv(2,0) = Jv(2,0) + (-dv) * Quad(1.d0,1.d0,1.d0,xc,xr) * ChiL_x(xc,xr) * ChiL_x(xc,xr)

     !     dL_v/dvr
     Jv(2,1) = Jv(2,1) + (-dv) * Quad(1.d0,1.d0,1.d0,xc,xr) * ChiL_x(xc,xr) * ChiR_x(xc,xr)

     !    L_w =  -dw * w_x * ChiL_x
     !     dL_w/dwc
     Jw(3,0) = Jw(3,0) + (-dw) * Quad(1.d0,1.d0,1.d0,xc,xr) * ChiL_x(xc,xr) * ChiL_x(xc,xr)

     !     dL_w/dwr
     Jw(3,1) = Jw(3,1) + (-dw) * Quad(1.d0,1.d0,1.d0,xc,xr) * ChiL_x(xc,xr) * ChiR_x(xc,xr)


     !    compute reaction Jacobian components

     !    R_u = (a - (w+1.d0)*u + v*u*u) 
     !     dR_u/duc
     df1 = (-(w1+1.d0) + 2.d0*v1*u1) * ChiL1 * ChiL1
     df2 = (-(w2+1.d0) + 2.d0*v2*u2) * ChiL2 * ChiL2
     df3 = (-(w3+1.d0) + 2.d0*v3*u3) * ChiL3 * ChiL3
     Ju(1,0) = Ju(1,0) + dQdf1*df1 + dQdf2*df2 + dQdf3*df3

     !     dR_u/dur 
     df1 = (-(w1+1.d0) + 2.d0*v1*u1) * ChiL1 * ChiR1
     df2 = (-(w2+1.d0) + 2.d0*v2*u2) * ChiL2 * ChiR2
     df3 = (-(w3+1.d0) + 2.d0*v3*u3) * ChiL3 * ChiR3
     Ju(1,1) = Ju(1,1) + dQdf1*df1 + dQdf2*df2 + dQdf3*df3

     !     dR_u/dvc
     df1 = (u1*u1) * ChiL1 * ChiL1
     df2 = (u2*u2) * ChiL2 * ChiL2
     df3 = (u3*u3) * ChiL3 * ChiL3
     Ju(2,0) = Ju(2,0) + dQdf1*df1 + dQdf2*df2 + dQdf3*df3

     !     dR_u/dvr 
     df1 = (u1*u1) * ChiL1 * ChiR1
     df2 = (u2*u2) * ChiL2 * ChiR2
     df3 = (u3*u3) * ChiL3 * ChiR3
     Ju(2,1) = Ju(2,1) + dQdf1*df1 + dQdf2*df2 + dQdf3*df3

     !     dR_u/dwc
     df1 = (-u1) * ChiL1 * ChiL1
     df2 = (-u2) * ChiL2 * ChiL2
     df3 = (-u3) * ChiL3 * ChiL3
     Ju(3,0) = Ju(3,0) + dQdf1*df1 + dQdf2*df2 + dQdf3*df3

     !     dR_u/dwr 
     df1 = (-u1) * ChiL1 * ChiR1
     df2 = (-u2) * ChiL2 * ChiR2
     df3 = (-u3) * ChiL3 * ChiR3
     Ju(3,1) = Ju(3,1) + dQdf1*df1 + dQdf2*df2 + dQdf3*df3


     !    R_v = (w*u - v*u*u) 
     !     dR_v/duc
     df1 = (w1 - 2.d0*v1*u1) * ChiL1 * ChiL1
     df2 = (w2 - 2.d0*v2*u2) * ChiL2 * ChiL2
     df3 = (w3 - 2.d0*v3*u3) * ChiL3 * ChiL3
     Jv(1,0) = Jv(1,0) + dQdf1*df1 + dQdf2*df2 + dQdf3*df3

     !     dR_v/dur 
     df1 = (w1 - 2.d0*v1*u1) * ChiL1 * ChiR1
     df2 = (w2 - 2.d0*v2*u2) * ChiL2 * ChiR2
     df3 = (w3 - 2.d0*v3*u3) * ChiL3 * ChiR3
     Jv(1,1) = Jv(1,1) + dQdf1*df1 + dQdf2*df2 + dQdf3*df3

     !     dR_v/dvc
     df1 = (-u1*u1) * ChiL1 * ChiL1
     df2 = (-u2*u2) * ChiL2 * ChiL2
     df3 = (-u3*u3) * ChiL3 * ChiL3
     Jv(2,0) = Jv(2,0) + dQdf1*df1 + dQdf2*df2 + dQdf3*df3

     !     dR_v/dvr 
     df1 = (-u1*u1) * ChiL1 * ChiR1
     df2 = (-u2*u2) * ChiL2 * ChiR2
     df3 = (-u3*u3) * ChiL3 * ChiR3
     Jv(2,1) = Jv(2,1) + dQdf1*df1 + dQdf2*df2 + dQdf3*df3

     !     dR_v/dwc
     df1 = (u1) * ChiL1 * ChiL1
     df2 = (u2) * ChiL2 * ChiL2
     df3 = (u3) * ChiL3 * ChiL3
     Jv(3,0) = Jv(3,0) + dQdf1*df1 + dQdf2*df2 + dQdf3*df3

     !     dR_v/dwr 
     df1 = (u1) * ChiL1 * ChiR1
     df2 = (u2) * ChiL2 * ChiR2
     df3 = (u3) * ChiL3 * ChiR3
     Jv(3,1) = Jv(3,1) + dQdf1*df1 + dQdf2*df2 + dQdf3*df3


     !    R_w = ((b-w)/ep - w*u) 
     !     dR_w/duc
     df1 = (-w1) * ChiL1 * ChiL1
     df2 = (-w2) * ChiL2 * ChiL2
     df3 = (-w3) * ChiL3 * ChiL3
     Jw(1,0) = Jw(1,0) + dQdf1*df1 + dQdf2*df2 + dQdf3*df3

     !     dR_w/dur 
     df1 = (-w1) * ChiL1 * ChiR1
     df2 = (-w2) * ChiL2 * ChiR2
     df3 = (-w3) * ChiL3 * ChiR3
     Jw(1,1) = Jw(1,1) + dQdf1*df1 + dQdf2*df2 + dQdf3*df3

     !     dR_w/dwc
     df1 = (-1.d0/ep - u1) * ChiL1 * ChiL1
     df2 = (-1.d0/ep - u2) * ChiL2 * ChiL2
     df3 = (-1.d0/ep - u3) * ChiL3 * ChiL3
     Jw(3,0) = Jw(3,0) + dQdf1*df1 + dQdf2*df2 + dQdf3*df3

     !     dR_w/dwr 
     df1 = (-1.d0/ep - u1) * ChiL1 * ChiR1
     df2 = (-1.d0/ep - u2) * ChiL2 * ChiR2
     df3 = (-1.d0/ep - u3) * ChiL3 * ChiR3
     Jw(3,1) = Jw(3,1) + dQdf1*df1 + dQdf2*df2 + dQdf3*df3


     ! insert Jacobian entries into CSR matrix structure

     !   Ju row
     Jrowptrs(idx(ix,1)+1) = nz

     Jdata(nz+1:nz+3) = (/ Ju(1,-1), Ju(2,-1), Ju(3,-1) /)
     Jcolvals(nz+1:nz+3) = (/ idx(ix-1,1), idx(ix-1,2), idx(ix-1,3) /)
     nz = nz+3

     Jdata(nz+1:nz+3) = (/ Ju(1,0), Ju(2,0), Ju(3,0) /)
     Jcolvals(nz+1:nz+3) = (/ idx(ix,1), idx(ix,2), idx(ix,3) /)
     nz = nz+3

     Jdata(nz+1:nz+3) = (/ Ju(1,1), Ju(2,1), Ju(3,1) /)
     Jcolvals(nz+1:nz+3) = (/ idx(ix+1,1), idx(ix+1,2), idx(ix+1,3) /)
     nz = nz+3

     !   Jv row
     Jrowptrs(idx(ix,2)+1) = nz

     Jdata(nz+1:nz+3) = (/ Jv(1,-1), Jv(2,-1), Jv(3,-1) /)
     Jcolvals(nz+1:nz+3) = (/ idx(ix-1,1), idx(ix-1,2), idx(ix-1,3) /)
     nz = nz+3

     Jdata(nz+1:nz+3) = (/ Jv(1,0), Jv(2,0), Jv(3,0) /)
     Jcolvals(nz+1:nz+3) = (/ idx(ix,1), idx(ix,2), idx(ix,3) /)
     nz = nz+3

     Jdata(nz+1:nz+3) = (/ Jv(1,1), Jv(2,1), Jv(3,1) /)
     Jcolvals(nz+1:nz+3) = (/ idx(ix+1,1), idx(ix+1,2), idx(ix+1,3) /)
     nz = nz+3

     !   Jw row
     Jrowptrs(idx(ix,3)+1) = nz

     Jdata(nz+1:nz+3) = (/ Jw(1,-1), Jw(2,-1), Jw(3,-1) /)
     Jcolvals(nz+1:nz+3) = (/ idx(ix-1,1), idx(ix-1,2), idx(ix-1,3) /)
     nz = nz+3

     Jdata(nz+1:nz+3) = (/ Jw(1,0), Jw(2,0), Jw(3,0) /)
     Jcolvals(nz+1:nz+3) = (/ idx(ix,1), idx(ix,2), idx(ix,3) /)
     nz = nz+3

     Jdata(nz+1:nz+3) = (/ Jw(1,1), Jw(2,1), Jw(3,1) /)
     Jcolvals(nz+1:nz+3) = (/ idx(ix+1,1), idx(ix+1,2), idx(ix+1,3) /)
     nz = nz+3


  enddo

  ! Dirichlet boundary at right
  Jrowptrs(idx(Nint,1)+1) = nz
  Jrowptrs(idx(Nint,2)+1) = nz
  Jrowptrs(idx(Nint,3)+1) = nz
 
  ! signal end of data in CSR matrix
  Jrowptrs(idx(Nint,3)+2) = nz

  ier = 0
  return

end subroutine farkspjac
!-----------------------------------------------------------------



subroutine farkspmass(t, neq, nnz, Mdata, Mcolvals, Mrowptrs, &
                      ipar, rpar, wk1, wk2, wk3, ier)
!-----------------------------------------------------------------
! Mass matrix computation routine
!-----------------------------------------------------------------
  use UserData
  use Quadrature
  use FEM

  ! Declarations
  implicit none

  ! Arguments
  real(kind=REALTYPE), intent(in)  :: t, rpar(1)
  real(kind=REALTYPE), intent(in), dimension(3,N) :: wk1, wk2, wk3
  real(kind=REALTYPE), intent(out) :: Mdata(nnz)
  integer*8, intent(in) :: ipar(1), neq, nnz
  integer(kind=SUNINDEXTYPE), intent(out) :: Mcolvals(nnz)
  integer(kind=SUNINDEXTYPE), intent(out) :: Mrowptrs(neq+1)
  integer,  intent(out) :: ier

  ! Local data
  integer :: ix, nz, Nint
  real*8  :: xl, xc, xr, Ml, Mc, Mr, ChiL1, ChiL2, ChiL3, ChiR1, ChiR2, ChiR3
  logical :: left, right

  ! check that vector/matrix dimensions match up
  if ((3*N /= neq) .or. (nnz /= 15*neq)) then
     ier = 1
     return
  endif

  ! set integer*4 version of N for call to idx()
  Nint = N
  
  ! clear out Jacobian matrix data
  Mdata = 0.d0
  nz = 0

  ! iterate through nodes, filling in matrix by rows
  do ix=1,N

     ! set booleans to determine whether intervals exist on the left/right */
     left  = .true.
     right = .true.
     if (ix==1)  left  = .false.
     if (ix==N)  right = .false.

     ! set nodal value shortcuts (interval index aligns with left node)
     if (left) then
        xl = x(ix-1)
     endif
     xc = x(ix)
     if (right) then
        xr = x(ix+1)
     endif

     ! compute entries of all mass matrix rows at node ix
     Ml = 0.d0
     Mc = 0.d0
     Mr = 0.d0

     ! first compute dependence on values to left and center
     if (left) then

        ChiL1 = ChiL(xl, xc, X1(xl,xc))
        ChiL2 = ChiL(xl, xc, X2(xl,xc))
        ChiL3 = ChiL(xl, xc, X3(xl,xc))
        ChiR1 = ChiR(xl, xc, X1(xl,xc))
        ChiR2 = ChiR(xl, xc, X2(xl,xc))
        ChiR3 = ChiR(xl, xc, X3(xl,xc))

        Ml = Ml + Quad(ChiL1*ChiR1, ChiL2*ChiR2, ChiL3*ChiR3, xl, xc)
        Mc = Mc + Quad(ChiR1*ChiR1, ChiR2*ChiR2, ChiR3*ChiR3, xl, xc)

     endif

     ! second compute dependence on values to center and right
     if (right) then

        ChiL1 = ChiL(xc, xr, X1(xc,xr))
        ChiL2 = ChiL(xc, xr, X2(xc,xr))
        ChiL3 = ChiL(xc, xr, X3(xc,xr))
        ChiR1 = ChiR(xc, xr, X1(xc,xr))
        ChiR2 = ChiR(xc, xr, X2(xc,xr))
        ChiR3 = ChiR(xc, xr, X3(xc,xr))

        Mc = Mc + Quad(ChiL1*ChiL1, ChiL2*ChiL2, ChiL3*ChiL3, xc, xr)
        Mr = Mr + Quad(ChiL1*ChiR1, ChiL2*ChiR2, ChiL3*ChiR3, xc, xr)

     endif


     ! insert mass matrix entries into CSR matrix structure

     !   u row
     Mrowptrs(idx(ix,1)+1) = nz
     if (left) then
        nz = nz+1
        Mdata(nz) = Ml
        Mcolvals(nz) = idx(ix-1,1)
     endif
     nz = nz+1
     Mdata(nz) = Mc
     Mcolvals(nz) = idx(ix,1)
     if (right) then
        nz = nz+1
        Mdata(nz) = Mr
        Mcolvals(nz) = idx(ix+1,1)
     endif

     !   v row
     Mrowptrs(idx(ix,2)+1) = nz
     if (left) then
        nz = nz+1
        Mdata(nz) = Ml
        Mcolvals(nz) = idx(ix-1,2)
     endif
     nz = nz+1
     Mdata(nz) = Mc
     Mcolvals(nz) = idx(ix,2)
     if (right) then
        nz = nz+1
        Mdata(nz) = Mr
        Mcolvals(nz) = idx(ix+1,2)
     endif

     !   w row
     Mrowptrs(idx(ix,3)+1) = nz
     if (left) then
        nz = nz+1
        Mdata(nz) = Ml
        Mcolvals(nz) = idx(ix-1,3)
     endif
     nz = nz+1
     Mdata(nz) = Mc
     Mcolvals(nz) = idx(ix,3)
     if (right) then
        nz = nz+1
        Mdata(nz) = Mr
        Mcolvals(nz) = idx(ix+1,3)
     endif

  enddo

  ! signal end of data in CSR matrix
  Mrowptrs(idx(Nint,3)+2) = nz

  ier = 0
  return

end subroutine farkspmass
!-----------------------------------------------------------------
