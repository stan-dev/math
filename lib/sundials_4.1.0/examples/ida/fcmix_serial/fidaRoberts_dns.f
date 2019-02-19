C     ----------------------------------------------------------------
C     Programmer(s): Daniel R. Reynolds @ SMU
C                    Alan C. Hindmarsh and Radu Serban @ LLNL      
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
c     This simple example problem for FIDA, due to Robertson, is from 
c     chemical kinetics, and consists of the following three equations:
c
c          dy1/dt = -.04*y1 + 1.e4*y2*y3
c          dy2/dt =  .04*y1 - 1.e4*y2*y3 - 3.e7*y2**2
c             0   = y1 + y2 + y3 - 1
c
c     on the interval from t = 0.0 to t = 4.e10, with initial
c     conditions: y1 = 1, y2 = y3 = 0.
c
c     While integrating the system, we also employ the rootfinding feature
c     to find the points at which y1 = 1.e-4 or at which y3 = 0.01.
c
c     The problem is solved using a dense linear solver, with a
c     user-supplied Jacobian. Output is printed at
c     t = .4, 4, 40, ..., 4e10.
c     ----------------------------------------------------------------
c
      program fidaRoberts_dns
c
      implicit none
c
c The following declaration specification should match C type long int.
      integer*8 neq, iout(25), ipar
      integer ier, ierroot, info(2)
      double precision rout(10), rpar
      integer iatol, nout, jout, itask
      integer nst, kused, hused, i
      double precision t0, t1, rtol, tout, tret
      double precision y(3), yp(3), atol(3)
c
      data nst/3/, kused/9/, hused/2/
c
c Initialize variables
c
      neq = 3
      nout = 12
      rtol = 1.0d-4
      t0 = 0.0d0
      t1 = 0.4d0
      iatol = 2
      itask = 1
c
      y(1) = 1.0d0
      y(2) = 0.0d0
      y(3) = 0.0d0
c
      yp(1) = -0.04d0
      yp(2) = 0.04d0
      yp(3) = 0.0d0
c
      atol(1) = 1.0d-6
      atol(2) = 1.0d-10
      atol(3) = 1.0d-6
c
c Initialize IDA vector environment
c
      call fnvinits(2, neq, ier)
      if (ier .ne. 0) then
         write(6,10) ier
 10      format(///' SUNDIALS_ERROR: FNVINITS returned IER = ', i5)
         stop
      endif
c
      call fidamalloc(t0, y, yp, iatol, rtol, atol, 
     &                iout, rout, ipar, rpar, ier)
      if (ier .ne. 0) then
         write(6,20) ier
 20      format(///' SUNDIALS_ERROR: FIDAMALLOC returned IER = ', i5)
         stop
      endif
c
c Initialize rootfinding problem

      call fidarootinit(2, ier)
      if (ier .ne. 0) then
         write(6,25) ier
 25      format(///' SUNDIALS_ERROR: FIDAROOTINIT returned IER = ', i5)
         call fidafree
         stop
      endif
c
c Attach dense matrix and linear solver
c
      call fsundensematinit(2, neq, neq, ier)
      if (ier .ne. 0) then
         write(6,30) ier
 30      format(///' SUNDIALS_ERROR: FSUNDENSEMATINIT IER = ', i5)
         call fidafree
         stop
      endif
      call fsundenselinsolinit(2, ier)
      if (ier .ne. 0) then
         write(6,33) ier
 33      format(///' SUNDIALS_ERROR: FSUNDENSELINSOLINIT IER = ', i5)
         call fidafree
         stop
      endif
      call fidalsinit(ier)
      if (ier .ne. 0) then
         write(6,35) ier
 35      format(///' SUNDIALS_ERROR: FIDALSINIT returned IER = ', i5)
         call fidafree
         stop
      endif
      call fidadensesetjac(1, ier)
      if (ier .ne. 0) then
         write(6,37) ier
 37      format(///' SUNDIALS_ERROR: FIDADENSESETJAC IER = ', i5)
         call fidafree
         stop
      endif
c
c Print header
c
      call prntintro(rtol, atol, y)
c
      tout = t1
c
c
      jout = 1
      do while(jout .le. nout)
c
        call fidasolve(tout, tret, y, yp, itask, ier)
c
        write(6,40) tret, (y(i), i = 1,3), iout(nst), iout(kused),
     &              rout(hused)
 40     format(e10.4, 3(1x,e12.4), i5, i3, e12.4)
c
        if (ier .lt. 0) then
           write(6,50) ier, iout(15)
 50        format(///' SUNDIALS_ERROR: FIDASOLVE returned IER = ',i5,/,
     1            '                 Linear Solver returned IER = ',i5)
           call fidarootfree
           call fidafree
           stop
        endif
c
        if (ier .eq. 2) then
          call fidarootinfo(2, info, ierroot)
          if (ierroot .lt. 0) then
            write(6,55) ierroot
 55         format(///' SUNDIALS_ERROR: FIDAROOTINFO returned IER = ',
     1             i5)
            call fidarootfree
            call fidafree
            stop
          endif
          write(6,60) (info(i), i = 1,2)
 60       format(5x, 'Above is a root, INFO() = ', 2i3)
        endif                   
c
        if (ier .eq. 0) then
           tout = tout * 10.0d0
           jout = jout + 1
        endif
c
      ENDDO
c
c Print final statistics
c
      call prntstats(iout)
c
c Free IDA memory
c
      call fidarootfree
      call fidafree
c
      stop
      end
c
c ==========
c
      subroutine fidaresfun(tres, y, yp, res, ipar, rpar, reserr)
c
      implicit none
c
c The following declaration specification should match C type long int.
      integer*8 ipar(*)
      integer reserr
      double precision tres, rpar(*)
      double precision y(*), yp(*), res(*)
c
      res(1) = -0.04d0*y(1)+1.0d4*y(2)*y(3)
      res(2) = -res(1)-3.0d7*y(2)*y(2)-yp(2)
      res(1) = res(1)-yp(1)
      res(3) = y(1)+y(2)+y(3)-1.0d0
c
      reserr = 0
c
      return
      end
c
c ==========
c
      subroutine fidadjac(neq, t, y, yp, r, jac, cj, ewt, h,
     1                    ipar, rpar, wk1, wk2, wk3, djacerr)
c
      implicit none
c
c The following declaration specification should match C type long int.
      integer*8 neq, ipar(*)
      integer djacerr
      double precision t, h, cj, rpar(*)
      double precision y(*), yp(*), r(*), ewt(*), jac(neq,neq)
      double precision wk1(*), wk2(*), wk3(*)
c
      jac(1,1) = -0.04d0-cj
      jac(2,1) = 0.04d0
      jac(3,1) = 1.0d0
      jac(1,2) = 1.0d4*y(3)
      jac(2,2) = -1.0d4*y(3)-6.0d7*y(2)-cj
      jac(3,2) = 1.0d0
      jac(1,3) = 1.0d4*y(2)
      jac(2,3) = -1.0d4*y(2)
      jac(3,3) = 1.0d0
c
      djacerr = 0

      return
      end
c
c ==========
c
      subroutine fidarootfn(t, y, yp, g, ipar, rpar, ier)
c Fortran routine for rootfinding
      implicit none
c
c The following declaration specification should match C type long int.
      integer*8 ipar(*)
      integer ier
      double precision t, y(*), yp(*), g(*), rpar(*)
c
      g(1) = y(1) - 1.0d-4
      g(2) = y(3) - 1.0d-2

      ier = 0

      return
      end
c
c ==========
c
      subroutine prntintro(rtol, atol, y)
c
      implicit none
c
      integer i
      double precision rtol, atol(*), y(*)
c
      write(6,60) rtol, (atol(i), i = 1,3), (y(i), i = 1,3)
 60   format(/'fidaRoberts_dns: Robertson kinetics DAE serial example',
     &       'problem for IDA', /,'          Three equation chemical',
     &       'kinetics problem.', //,
     &       'Tolerance parameters:  rtol = ', e8.2,
     &       '   atol = ', 3(1x,e8.2), /,
     &       'Initial conditions y0 = (', 3(1x,e8.2), ')', //,
     &       '  t            y1           y2           y3        nst',
     &       '  k    h')
c
      return
      end
c
c ==========
c
      subroutine prntstats(iout)
c
      implicit none
c
c The following declaration specification should match C type long int.
      integer*8 iout(25)
      integer nst, reseval, jaceval, nni, ncf, netf, nge
c
      data nst/3/, reseval/4/, jaceval/17/, nni/7/, netf/5/,
     &     ncf/6/, nge/12/
c
      write(6,70) iout(nst), iout(reseval), iout(jaceval),
     &             iout(nni), iout(netf), iout(ncf), iout(nge)
 70   format(/'Final Run Statistics:', //,
     &         'Number of steps                    = ', i3, /,
     &         'Number of residual evaluations     = ', i3, /,
     &         'Number of Jacobian evaluations     = ', i3, /,
     &         'Number of nonlinear iterations     = ', i3, /,
     &         'Number of error test failures      = ', i3, /,
     &         'Number of nonlinear conv. failures = ', i3, /,
     &         'Number of root function evals.     = ', i3)
c
      return
      end
