C     ----------------------------------------------------------------
C     Programmer(s): Daniel R. Reynolds @ SMU
C                 Radu Serban and Alan C. Hindmarsh @ LLNL      
C     ----------------------------------------------------------------
C     SUNDIALS Copyright Start
C     Copyright (c) 2002-2019, Lawrence Livermore National Security
C     and Southern Methodist University.
C     All rights reserved.
C
C     See the top-level LICENSE and NOTICE files for details.
C
C     SPDX-License-Identifier: BSD-3-Clause
C     SUNDIALS Copyright End
c     ----------------------------------------------------------------
c     Example problem for FIDA: 2D heat equation, parallel, GMRES,
c     IDABBDPRE.
c
c     This example solves a discretized 2D heat equation problem.
c     This version uses the Krylov solver SUNSPGMR and BBD
c     preconditioning.
c
c     The DAE system solved is a spatial discretization of the PDE
c              du/dt = d^2u/dx^2 + d^2u/dy^2
c     on the unit square. The boundary condition is u = 0 on all edges.
c     Initial conditions are given by u = 16 x (1 - x) y (1 - y). The
c     PDE is treated with central differences on a uniform MX x MY
c     grid. The values of u at the interior points satisfy ODEs, and
c     equations u = 0 at the boundaries are appended, to form a DAE
c     system of size N = MX * MY. Here MX = MY = 10.
c
c     The system is actually implemented on submeshes, processor by
c     processor, with an MXSUB by MYSUB mesh on each of NPEX * NPEY
c     processors.
c
c     The system is solved with FIDA using the Krylov linear solver
c     SUNSPGMR in conjunction with the preconditioner module IDABBDPRE.
c     The preconditioner uses a tridiagonal approximation
c     (half-bandwidths = 1). The constraints u >= 0 are posed for all
c     components. Local error testing on the boundary values is
c     suppressed. Output is taken at t = 0, .01, .02, .04, ..., 10.24.
c     ----------------------------------------------------------------
c
      program fidaHeat2D_kry_bbd_p
c
      include "mpif.h"
c
c global variables
c
c The following declaration specification should match C type long int.
      integer*8 nlocal, neq
      integer npex, npey, mxsub, mysub, mx, my, ixsub, jysub
      integer thispe, mxsubg, mysubg
      parameter (mxsubg = 5, mysubg = 5)
      double precision dx, dy, coeffx, coeffy, coeffxy
      double precision uext((mxsubg+2)*(mysubg+2))
c
c local variables
c
c The following declaration specification should match C type long int.
      integer*8 mudq, mldq, mukeep, mlkeep, iout(25), ipar(1)
      double precision rout(10), rpar(1)
      integer nout, ier
      parameter (nout = 11)
      integer npes, inopt, maxl, maxrs, itask, iatol
      integer nlocalg
      parameter (nlocalg = mxsubg*mysubg)
      double precision t0, t1, tout, tret, dqrely
      double precision atol, rtol
      double precision constr(nlocalg), uu(nlocalg), up(nlocalg)
      double precision res(nlocalg), id(nlocalg)
c
      data atol/1.0d-3/, rtol/0.0d0/
c
      common /pcom/ dx, dy, coeffx, coeffy, coeffxy, uext,
     &              nlocal, neq, mx, my, mxsub, mysub, npex, npey,
     &              ixsub, jysub, thispe
c
c Initialize variables
c
      npex = 2
      npey = 2
      mxsub = 5
      mysub = 5
      mx = npex*mxsub
      my = npey*mysub
      neq = mx*my
      nlocal = mxsub*mysub
      inopt = 1
      t0 = 0.0d0
      t1 = 0.01d0
      mudq = mxsub
      mldq = mxsub
      mukeep = 1
      mlkeep = 1
      dqrely = 0.0d0
      maxl = 0
      maxrs = 0
      itask = 1
      iatol = 1
c
c Initialize MPI environment
c
      call mpi_init(ier)
      if (ier .ne. 0) then
         write(*,2) ier
 2       format(///' MPI_ERROR: MPI_INIT returned IER = ', i5)
         stop
      endif
c
      call mpi_comm_size(mpi_comm_world, npes, ier)
      if (ier .ne. 0) then
         write(*,3) ier
 3       format(///' MPI_ERROR: MPI_COMM_SIZE returned IER = ', i5)
         call mpi_abort(mpi_comm_world, 1, ier)
         stop
      endif
c
      call mpi_comm_rank(mpi_comm_world, thispe, ier)
      if (ier .ne. 0) then
         write(*,4) ier
 4       format(///' MPI_ERROR: MPI_COMM_RANK returned IER = ', i5)
         call mpi_abort(mpi_comm_world, 1, ier)
         stop
      endif
c
      if (npes .ne. npex*npey) then
         if (thispe .eq. 0) then
            write(*,5) npes, npex*npey
 5          format(///' MPI_ERROR: npes = ', i5, ' is not equal to ',
     &                'NPEX*NPEY = ', i5)
            call mpi_finalize(ier)
            stop
         endif
      endif
c
      call fnvinitp(mpi_comm_world, 2, nlocal, neq, ier)
      if (ier .ne. 0) then
         write(*,6) ier
 6       format(///' SUNDIALS_ERROR: FNVINITP returned IER = ', i5)
         call mpi_finalize(ier)
         stop
      endif
c
      jysub = int(thispe/npex)
      ixsub = thispe-jysub*npex
c
c Initialize problem data
c
      call setinitprofile(uu, up, id, res, constr, ipar, rpar)
c
c Initialize IDA environment
c
      call fidamalloc(t0, uu, up, iatol, rtol, atol, 
     &     iout, rout, ipar, rpar, ier)
      if (ier .ne. 0) then
         write(*,7) ier
 7       format(///' SUNDIALS_ERROR: FIDAMALLOC returned IER = ', i5)
         call mpi_abort(mpi_comm_world, 1, ier)
         stop
      endif
c
c     Set optional inputs
c
      call fidasetiin('SUPPRESS_ALG', 1, ier)
      call fidasetvin('ID_VEC', id, ier)
      call fidasetvin('CONSTR_VEC', constr, ier)
c
c Initialize and attach SUNSPGMR solver
c
      call fsunspgmrinit(2, 1, maxl, ier)
      if (ier .ne. 0) then
         write(*,8) ier
 8       format(///' SUNDIALS_ERROR: FSUNSPGMRINIT returned IER = ', i5)
         call mpi_abort(mpi_comm_world, 1, ier)
         stop
      endif
      call fsunspgmrsetmaxrs(2, 5, ier)
      if (ier .ne. 0) then
         write(*,9) ier
 9       format(///' SUNDIALS_ERROR: FSUNSPGMRSETMAXRS IER = ', i5)
         call mpi_abort(mpi_comm_world, 1, ier)
         stop
      endif
      call fidalsinit(ier)
      if (ier .ne. 0) then
         write(*,10) ier
 10       format(///' SUNDIALS_ERROR: FIDALSINIT IER = ', i5)
         call mpi_abort(mpi_comm_world, 1, ier)
         stop
      endif
c
c Initialize and attach BBD preconditioner module
c
      call fidabbdinit(nlocal, mudq, mldq, mukeep, mlkeep, dqrely, ier)
      if (ier .ne. 0) then
         write(*,11) ier
 11       format(///' SUNDIALS_ERROR: FIDABBDINIT returned IER = ', i5)
         call mpi_abort(mpi_comm_world, 1, ier)
         stop
      endif
c
c Print header
c
      if (thispe .eq. 0) then
         call prntintro(rtol, atol)
         call prntcase(1, mudq, mukeep)
      endif
c
      tout = t1
      do 12 jout = 1, nout
c
         call fidasolve(tout, tret, uu, up, itask, ier)
c
         call prntoutput(tret, uu, iout, rout)
c
         if (ier .ne. 0) then
            write(*,13) ier
 13         format(///' SUNDIALS_ERROR: FIDASOLVE returned IER = ', i5)
            call fidafree
            stop
         endif
c
         tout = tout*2.0d0
c
 12   continue
c
c Print statistics
c
      if (thispe .eq. 0) then
         call prntfinalstats(iout)
      endif
c
c Reinitialize variables and data for second problem
c
      mudq = 1
      mldq = 1
c
      call setinitprofile(uu, up, id, res, constr, ipar, rpar)
c
      call fidareinit(t0, uu, up, iatol, rtol, atol, ier)
      if (ier .ne. 0) then
         write(*,43) ier
 43      format(///' SUNDIALS_ERROR: FIDAREINIT returned IER = ', i5)
      endif
c
      call fidabbdreinit(nlocal, mudq, mldq, dqrely, ier)
      if (ier .ne. 0) then
         write(*,44) ier
 44      format(///' SUNDIALS_ERROR: FIDABBDREINIT returned IER = ', i5)
         call fidafree
         stop
      endif
c
c Print header
c
      if (thispe .eq. 0) then
         call prntcase(2, mudq, mukeep)
      endif
c
      tout = t1
      do 14 jout = 1, nout
c
         call fidasolve(tout, tret, uu, up, itask, ier)
c
         call prntoutput(tret, uu, iout, rout)
c
         if (ier .ne. 0) then
            write(*,15) ier
 15         format(///' SUNDIALS_ERROR: FIDASOLVE returned IER = ', i5)
            call fidafree
            stop
         endif
c
         tout = tout*2.0d0
c
 14   continue
c
c Print statistics
c
      if (thispe .eq. 0) then
         call prntfinalstats(iout)
      endif
c
c Free memory
c
      call fidafree
c
      call mpi_finalize(ier)
c
      stop
      end
c
c ==========
c
      subroutine setinitprofile(uu, up, id, res, constr, ipar, rpar)
c
c global variables
c
c The following declaration specification should match C type long int.
      integer*8 nlocal, neq
      integer npex, npey, mxsub, mysub, mx, my, ixsub, jysub
      integer thispe, mxsubg, mysubg
      parameter (mxsubg = 5, mysubg = 5)
      double precision dx, dy, coeffx, coeffy, coeffxy, rpar(*)
      double precision uext((mxsubg+2)*(mysubg+2))
c
c local variables
c
c The following declaration specification should match C type long int.
      integer*8 ipar(*)
      integer i, iloc, j, jloc, offset, loc
      integer ixbegin, ixend, jybegin, jyend
      integer reserr
      double precision xfact, yfact
      double precision uu(*), up(*), id(*), res(*), constr(*)
c
      common /pcom/ dx, dy, coeffx, coeffy, coeffxy, uext,
     &              nlocal, neq, mx, my, mxsub, mysub, npex, npey,
     &              ixsub, jysub, thispe
c
c Initialize variables
c
      dx = 1.0d0/dble(mx-1)
      dy = 1.0d0/dble(my-1)
      coeffx = 1.0d0/(dx*dx)
      coeffy = 1.0d0/(dy*dy)
      coeffxy = 2.0d0/(dx*dx)+2.0d0/(dy*dy)
      ixbegin = mxsub*ixsub
      ixend = mxsub*(ixsub+1)-1
      jybegin = mysub*jysub
      jyend = mysub*(jysub+1)-1
c
      do 16 i = 1, nlocal
         id(i) = 1.0d0
 16   continue
c
      jloc = 0
      do 17 j = jybegin, jyend
         yfact = dy*dble(j)
         offset = jloc*mxsub
         iloc = 0
         do 18 i = ixbegin, ixend
            xfact = dx*dble(i)
            loc = offset+iloc
            uu(loc+1) = 16.0d0*xfact*(1.0d0-xfact)*yfact*(1.0d0-yfact)
            if (i .eq. 0 .or. i .eq. mx-1) then
               id(loc+1) = 0.0d0
            endif
            if (j .eq. 0 .or. j .eq. my-1) then
               id(loc+1) = 0.0d0
            endif
            iloc = iloc+1
 18      continue
         jloc = jloc+1
 17   continue
c
      do 19 i = 1, nlocal
         up(i) = 0.0d0
         constr(i) = 1.0d0
 19   continue
c
      call fidaresfun(0.0d0, uu, up, res, ipar, rpar, reserr)
c
      do 20 i = 1, nlocal
         up(i) = -1.0d0*res(i)
 20   continue
c
      return
      end
c
c ==========
c
      subroutine fidaresfun(tres, u, up, res, ipar, rpar, reserr)
c
c global variables
c
c The following declaration specification should match C type long int.
      integer*8 nlocal, neq
      integer npex, npey, mxsub, mysub, mx, my, ixsub, jysub
      integer thispe, mxsubg, mysubg
      parameter (mxsubg = 5, mysubg = 5)
      double precision dx, dy, coeffx, coeffy, coeffxy, rpar(*)
      double precision uext((mxsubg+2)*(mysubg+2))
c
c local variables
c
c The following declaration specification should match C type long int.
      integer*8 ipar(*)
      integer reserr
      double precision tres
      double precision u(*), up(*), res(*)
c
      common /pcom/ dx, dy, coeffx, coeffy, coeffxy, uext,
     &              nlocal, neq, mx, my, mxsub, mysub, npex, npey,
     &              ixsub, jysub, thispe
c
      call fidacommfn(nlocal, tres, u, up, ipar, rpar, reserr)
c
      call fidaglocfn(nlocal, tres, u, up, res, ipar, rpar, reserr)
c
      return
      end
c
c ==========
c
      subroutine fidacommfn(nloc, tres, u, up, ipar, rpar, reserr)
c
      include "mpif.h"
c
c global variables
c
c The following declaration specification should match C type long int.
      integer*8 nlocal, neq
      integer npex, npey, mxsub, mysub, mx, my, ixsub, jysub
      integer thispe, mxsubg, mysubg
      parameter (mxsubg = 5, mysubg = 5)
      double precision dx, dy, coeffx, coeffy, coeffxy, rpar(*)
      double precision uext((mxsubg+2)*(mysubg+2))
c
c local variables
c
c The following declaration specification should match C type long int.
      integer*8 nloc, ipar(*)
      integer reserr
      double precision tres, u(*), up(*)
c
      integer request(mpi_status_size)
      double precision buffer(2*mysub)
c
      common /pcom/ dx, dy, coeffx, coeffy, coeffxy, uext,
     &              nlocal, neq, mx, my, mxsub, mysub, npex, npey,
     &              ixsub, jysub, thispe
c
      call brecvpost(request, mxsub, mysub, buffer)
c
      call bsend(mxsub, mysub, u)
c
      call brecvwait(request, mxsub, buffer)
c
      return
      end
c
c ==========
c
      subroutine fidaglocfn(nloc, tres, u, up, res, ipar, rpar, reserr)
c
c global variables
c
c The following declaration specification should match C type long int.
      integer*8 nlocal, neq
      integer npex, npey, mxsub, mysub, mx, my, ixsub, jysub
      integer thispe, mxsubg, mysubg
      parameter (mxsubg = 5, mysubg = 5)
      double precision dx, dy, coeffx, coeffy, coeffxy, rpar(*)
      double precision uext((mxsubg+2)*(mysubg+2))
c
c local variables
c
c The following declaration specification should match C type long int.
      integer*8 nloc, ipar(*)
      integer reserr
      double precision tres, u(*), up(*), res(*)
c
      integer i, lx, ly, offsetu, offsetue, locu, locue
      integer ixbegin, ixend, jybegin, jyend, mxsub2
      double precision termx, termy, termctr
c
      common /pcom/ dx, dy, coeffx, coeffy, coeffxy, uext,
     &              nlocal, neq, mx, my, mxsub, mysub, npex, npey,
     &              ixsub, jysub, thispe
c
      mxsub2 = mxsub+2
c
      do 21 i = 1, nlocal
         res(i) = u(i)
 21   continue
c
      offsetu = 0
      offsetue = mxsub2+1
      do 22 ly = 0, mysub-1
         do 23 lx = 0, mxsub-1
            uext(offsetue+lx+1) = u(offsetu+lx+1)
 23      continue
         offsetu = offsetu+mxsub
         offsetue = offsetue+mxsub2
 22   continue
c
      ixbegin = 0
      ixend = mxsub-1
      jybegin = 0
      jyend = mysub-1
      if (ixsub .eq. 0) then
         ixbegin = ixbegin+1
      endif
      if (ixsub .eq. npex-1) then
         ixend = ixend-1
      endif
      if (jysub .eq. 0) then
         jybegin = jybegin+1
      endif
      if (jysub .eq. npey-1) then
         jyend = jyend-1
      endif
c
      do 24 ly = jybegin, jyend
         do 25 lx = ixbegin, ixend
            locu = lx+ly*mxsub
            locue = (lx+1)+(ly+1)*mxsub2
            termx = coeffx*(uext(locue)+uext(locue+2))
            termy = coeffy*(uext(locue-mxsub2+1)+uext(locue+mxsub2+1))
            termctr = coeffxy*uext(locue+1)
            res(locu+1) = up(locu+1)-(termx+termy-termctr)
 25      continue
 24   continue
c
      return
      end
c
c ==========
c
      subroutine bsend(dsizex, dsizey, uarray)
c
      include "mpif.h"
c
c global variables
c
c The following declaration specification should match C type long int.
      integer*8 nlocal, neq
      integer npex, npey, mxsub, mysub, mx, my, ixsub, jysub
      integer thispe, mxsubg, mysubg
      parameter (mxsubg = 5, mysubg = 5)
      double precision dx, dy, coeffx, coeffy, coeffxy
      double precision uext((mxsubg+2)*(mysubg+2))
c
c local variables
c
      integer dsizex, dsizey
      double precision uarray(*)
c
      integer ier, offsetu
      double precision bufleft(mysub), bufright(mysub)
c
      common /pcom/ dx, dy, coeffx, coeffy, coeffxy, uext,
     &              nlocal, neq, mx, my, mxsub, mysub, npex, npey,
     &              ixsub, jysub, thispe
c
      if (jysub .ne. 0) then
         call mpi_send(uarray(1), dsizex, mpi_double_precision,
     &                 thispe-npex, 0, mpi_comm_world, ier)
      endif
c
      if (jysub .ne. npey-1) then
         offsetu = (mysub-1)*dsizex
         call mpi_send(uarray(offsetu+1), dsizex, mpi_double_precision,
     &                 thispe+npex, 0, mpi_comm_world, ier)
      endif
c
      if (ixsub .ne. 0) then
         do 26 ly = 0, mysub-1
            offsetu = ly*dsizex
            bufleft(ly+1) = uarray(offsetu+1)
 26      continue
         call mpi_send(bufleft(1), dsizey, mpi_double_precision,
     &                 thispe-1, 0, mpi_comm_world, ier)
      endif
c
      if (ixsub .ne. npex-1) then
         do 27 ly = 0, mysub-1
            offsetu = ly*mxsub+(mxsub-1)
            bufright(ly+1) = uarray(offsetu+1)
 27      continue
         call mpi_send(bufright(1), dsizey, mpi_double_precision,
     &                 thispe+1, 0, mpi_comm_world, ier)
      endif
c
      return
      end
c
c ==========
c
      subroutine brecvpost(request, dsizex, dsizey, buffer)
c
      include "mpif.h"
c
c global variables
c
c The following declaration specification should match C type long int.
      integer*8 nlocal, neq
      integer npex, npey, mxsub, mysub, mx, my, ixsub, jysub
      integer thispe, mxsubg, mysubg
      parameter (mxsubg = 5, mysubg = 5)
      double precision dx, dy, coeffx, coeffy, coeffxy
      double precision uext((mxsubg+2)*(mysubg+2))
c
c local variables
c
      integer dsizex, dsizey
      integer request(*)
      double precision buffer(*)
      integer ier
      integer offsetue
c
      common /pcom/ dx, dy, coeffx, coeffy, coeffxy, uext,
     &              nlocal, neq, mx, my, mxsub, mysub, npex, npey,
     &              ixsub, jysub, thispe
c
      if (jysub .ne. 0) then
         call mpi_irecv(uext(2), dsizex, mpi_double_precision,
     &                  thispe-npex, 0, mpi_comm_world, request(1),
     &                  ier)
      endif
c
      if (jysub .ne. npey-1) then
         offsetue = (1+(mysub+1)*(mxsub+2))
         call mpi_irecv(uext(offsetue+1), dsizex, mpi_double_precision,
     &                  thispe+npex, 0, mpi_comm_world, request(2),
     &                  ier)
      endif
c
      if (ixsub .ne. 0) then
         call mpi_irecv(buffer(1), dsizey, mpi_double_precision,
     &                  thispe-1, 0, mpi_comm_world, request(3),
     &                  ier)
      endif
c
      if (ixsub .ne. npex-1) then
         call mpi_irecv(buffer(1+mysub), dsizey, mpi_double_precision,
     &                  thispe+1, 0, mpi_comm_world, request(4),
     &                  ier)
      endif
c
      return
      end
c
c ==========
c
      subroutine brecvwait(request, dsizex, buffer)
c
      include "mpif.h"
c
c global variables
c
c The following declaration specification should match C type long int.
      integer*8 nlocal, neq
      integer npex, npey, mxsub, mysub, mx, my, ixsub, jysub
      integer thispe, mxsubg, mysubg
      parameter (mxsubg = 5, mysubg = 5)
      double precision dx, dy, coeffx, coeffy, coeffxy
      double precision uext((mxsubg+2)*(mysubg+2))
c
c local variables
c
      integer request(*)
      integer dsizex
      double precision buffer(*)
      integer ly, dsizex2, offsetue
      integer ier, status(mpi_status_size)
c
      common /pcom/ dx, dy, coeffx, coeffy, coeffxy, uext,
     &              nlocal, neq, mx, my, mxsub, mysub, npex, npey,
     &              ixsub, jysub, thispe
c
      dsizex2 = dsizex+2
c
      if (jysub .ne. 0) then
         call mpi_wait(request(1), status, ier)
      endif
c
      if (jysub .ne. npey-1) then
         call mpi_wait(request(2), status, ier)
      endif
c
      if (ixsub .ne. 0) then
         call mpi_wait(request(3), status, ier)
         do 28 ly = 0, mysub-1
            offsetue = (ly+1)*dsizex2
            uext(offsetue+1) = buffer(ly+1)
 28      continue
      endif
c
      if (ixsub .ne. npex-1) then
         call mpi_wait(request(4), status, ier)
         do 29 ly = 0, mysub-1
            offsetue = (ly+2)*dsizex2-1
            uext(offsetue+1) = buffer(ly+mysub+1)
 29      continue
      endif
c
      return
      end
c
c ==========
c
      subroutine prntoutput(tret, u, iout, rout)
c
c global variables
c
c The following declaration specification should match C type long int.
      integer*8 nlocal, neq
      integer npex, npey, mxsub, mysub, mx, my, ixsub, jysub
      integer thispe, mxsubg, mysubg
      parameter (mxsubg = 5, mysubg = 5)
      double precision dx, dy, coeffx, coeffy, coeffxy
      double precision uext((mxsubg+2)*(mysubg+2))
c
c  local variables
c
c The following declaration specification should match C type long int.
      integer*8 iout(*), lenrwbbd, leniwbbd, ngebbd
      double precision tret, umax, u(*), rout(*)
c
      common /pcom/ dx, dy, coeffx, coeffy, coeffxy, uext,
     &              nlocal, neq, mx, my, mxsub, mysub, npex, npey,
     &              ixsub, jysub, thispe
c
      call maxnorm(u, umax)
c
      if (thispe .eq. 0) then
         call fidabbdopt(lenrwbbd, leniwbbd, ngebbd)
         write(*,30) tret, umax, iout(9), iout(3), iout(7),
     &               iout(22), iout(4), iout(16), ngebbd, rout(2),
     &               iout(20), iout(21)
 30      format(' ', e10.4, ' ', e13.5, '  ', i1, '  ', i2,
     &          '  ', i3, '  ', i3, '  ', i2,'+',i2, '  ',
     &          i3, '  ', e9.2, '  ', i2, '  ', i3)
      endif
c
      return
      end
c
c ==========
c
      subroutine maxnorm(u, unorm)
c
      include "mpif.h"
c
c global variables
c
c The following declaration specification should match C type long int.
      integer*8 nlocal, neq
      integer npex, npey, mxsub, mysub, mx, my, ixsub, jysub
      integer thispe, mxsubg, mysubg
      parameter (mxsubg = 5, mysubg = 5)
      double precision dx, dy, coeffx, coeffy, coeffxy
      double precision uext((mxsubg+2)*(mysubg+2))
c
c local variables
c
      integer i, ier
      double precision temp, unorm, u(*)
c
      common /pcom/ dx, dy, coeffx, coeffy, coeffxy, uext,
     &              nlocal, neq, mx, my, mxsub, mysub, npex, npey,
     &              ixsub, jysub, thispe
c
      temp = 0.0d0
c
      do 31 i = 1, nlocal
         temp = max(abs(u(i)), temp)
 31   continue
c
      call mpi_allreduce(temp, unorm, 1, mpi_double_precision,
     &                   mpi_max, mpi_comm_world, ier)
c
c      unorm = temp
c
      return
      end
c
c ==========
c
      subroutine prntintro(rtol, atol)
c
c global variables
c
c The following declaration specification should match C type long int.
      integer*8 nlocal, neq
      integer npex, npey, mxsub, mysub, mx, my, ixsub, jysub
      integer thispe, mxsubg, mysubg
      parameter (mxsubg = 5, mysubg = 5)
      double precision dx, dy, coeffx, coeffy, coeffxy
      double precision uext((mxsubg+2)*(mysubg+2))
c
c local variables
c
      double precision rtol, atol
c
      common /pcom/ dx, dy, coeffx, coeffy, coeffxy, uext,
     &              nlocal, neq, mx, my, mxsub, mysub, npex, npey,
     &              ixsub, jysub, thispe
c
      write(*,32) mx, my, neq, mxsub, mysub, npex, npey, rtol, atol
 32   format(/'fidaHeat2D_kry_bbd_p: Heat equation, parallel example',
     &     ' for FIDA', /, 16x,'Discretized heat equation',
     &     ' on 2D unit square.', /, 16x,'Zero boundary',
     &     ' conditions, polynomial conditions.', /,
     &     16x,'Mesh dimensions: ', i2, ' x ', i2,
     &     '         Total system size: ', i3, //,
     &     'Subgrid dimensions: ', i2, ' x ', i2,
     &     '           Processor array: ', i2, ' x ', i2, /,
     &     'Tolerance parameters: rtol = ', e8.2, '   atol = ',
     &     e8.2, /, 'Constraints set to force all solution',
     &     ' components >= 0.', /, 'SUPPRESSALG = SUNTRUE to remove',
     &     ' boundary components from the error test.', /,
     &     'Linear solver: SPGMR.    Preconditioner: BBDPRE - ',
     &     'Banded-block-diagonal.')
c     
      return
      end
c
c ==========
c
      subroutine prntcase(num, mudq, mukeep)
c
c local variables
c
c The following declaration specification should match C type long int.
      integer*8 mudq, mukeep
      integer num
c
      write(*,33) num, mudq, mukeep
 33   format(//, 'Case ', i2, /, '  Difference quotient half-',
     &     'bandwidths =', i2, /, '  Retained matrix half-bandwidths =',
     &     i2, //, 'Output Summary',/,'  umax = max-norm of solution',
     &     /,'  nre = nre + nreLS (total number of RES evals.)',
     &     //, '   time         umax       k nst  nni  nli   nre',
     &     '   nge      h     npe  nps', /,
     &     '-------------------------------------------------------',
     &     '-------------------')
c
      return
      end
c
c ==========
c
      subroutine prntfinalstats(iout)
c
c local variables
c
c The following declaration specification should match C type long int.
      integer*8 iout(*)
c
      write(*,34) iout(5), iout(6), iout(23)
 34   format(/, 'Error test failures            =', i3, /,
     &     'Nonlinear convergence failures =', i3, /,
     &     'Linear convergence failures    =', i3)
c
      return
      end
