C     --------------------------------------------------------------------
C     SUNDIALS Copyright Start
C     Copyright (c) 2002-2019, Lawrence Livermore National Security
C     and Southern Methodist University.
C     All rights reserved.
C
C     See the top-level LICENSE and NOTICE files for details.
C
C     SPDX-License-Identifier: BSD-3-Clause
C     SUNDIALS Copyright End
C     --------------------------------------------------------------------
C     FCVODE Example Problem: Advection-diffusion, band linear solver
C                             with banded user Jacobian.
C
C     The following is a simple example problem with a banded
C     Jacobian. The problem is the semi-discrete form of the
C     advection-diffusion equation in 2D:
C     du/dt = d^2 u / dx^2 + .5 du/dx + d^2 u / dy^2
C     on the rectangle 0 <= x <= 2, 0 <= y <= 1, and the time
C     interval 0 <= t <= 1. Homogeneous Dirichlet boundary conditions
C     are posed, and the initial condition is the following:
C     u(x,y,t=0) = x(2-x)y(1-y)exp(5xy) .
C     The PDE is discretized on a uniform MX+2 by MY+2 grid with
C     central differencing, and with boundary values eliminated,
C     leaving an ODE system of size NEQ = MX*MY.
C     This program solves this problem with CVODE, using the Fortran/C 
C     interface routine package. This solution uses the BDF method,
C     a user-supplied banded Jacobian routine, and scalar relative and
C     absolute tolerances. It prints results at t = .1, .2, ..., 1.0.
C     At the end of the run, various counters of interest are printed.
C     --------------------------------------------------------------------
C
      IMPLICIT NONE
C
      INTEGER*4 MX, MY, MXMY
      PARAMETER (MX=10, MY=5)
      PARAMETER (MXMY=MX*MY)
C     
      DOUBLE PRECISION XMAX, YMAX
      DATA XMAX/2.0D0/, YMAX/1.0D0/
C
      INTEGER*4 LNST, LNFE, LNSETUP, LNNI, LNCF, LNETF, LNJE
      INTEGER*4 IER, METH, IATOL, ITASK, JOUT
C The following declaration specification should match C type long int.
      INTEGER*8 NEQ, IOUT(25), IPAR(2), MU, ML
      DOUBLE PRECISION RTOL, ATOL, T0, T, TOUT, DTOUT, UNORM 
      DOUBLE PRECISION U(MXMY), ROUT(10), RPAR(5)
C
      DATA LNST/3/, LNFE/4/, LNETF/5/,  LNCF/6/, LNNI/7/, LNSETUP/8/, 
     1     LNJE/17/
C
      CALL INITBX(XMAX, YMAX, MX, MY, U, IPAR, RPAR)
C
      NEQ = MX*MY
      T0 = 0.0D0
      METH = 2
      IATOL = 1
      RTOL = 0.0D0
      ATOL = 1.0D-5
      MU = MY
      ML = MY
      DTOUT = 0.1D0
      ITASK = 1
C
      WRITE(6,10) NEQ
 10   FORMAT('Band example problem:'//
     1       ' Advection-diffusion, NEQ = ', I2//)
C
      CALL FNVINITS(1, NEQ, IER)
      IF (IER .NE. 0) THEN
        WRITE(6,20) IER
 20     FORMAT(///' SUNDIALS_ERROR: FNVINITS returned IER = ', I5)
        STOP
      ENDIF
C
C     initialize banded matrix module
      call FSUNBANDMATINIT(1, NEQ, MU, ML, MU+ML, IER)
      IF (IER .NE. 0) THEN
        WRITE(6,25) IER
 25     FORMAT(///' SUNDIALS_ERROR: FSUNBANDEMATINIT IER = ', I5)
        STOP
      ENDIF
C
C     initialize banded linear solver module
      call FSUNBANDLINSOLINIT(1, IER)
      IF (IER .NE. 0) THEN
        WRITE(6,28) IER
 28     FORMAT(///' SUNDIALS_ERROR: FSUNBANDLINSOLINIT IER = ', I5)
        STOP
      ENDIF
C
      CALL FCVMALLOC(T0, U, METH, IATOL, RTOL, ATOL,
     1               IOUT, ROUT, IPAR, RPAR, IER)
      IF (IER .NE. 0) THEN
        WRITE(6,30) IER
 30     FORMAT(///' SUNDIALS_ERROR: FCVMALLOC returned IER = ', I5)
        STOP
      ENDIF
C
C     attach matrix and linear solver modules to CVLs interface
      CALL FCVLSINIT(IER)
      IF (IER .NE. 0) THEN
        WRITE(6,40) IER
 40     FORMAT(///' SUNDIALS_ERROR: FCVLSINIT returned IER = ',I5)
        CALL FCVFREE
        STOP
      ENDIF
C
      CALL FCVBANDSETJAC(1, IER)
      IF (IER .NE. 0) THEN
        WRITE(6,45) IER
 45     FORMAT(///' SUNDIALS_ERROR: FCVDENSESETJAC returned IER = ',I5)
        CALL FCVFREE
        STOP
      ENDIF
C
      CALL MAXNORM(NEQ, U, UNORM)
      WRITE(6,48) T0, UNORM
 48   FORMAT(' At t = ', F6.2, '  max.norm(u) = ', E14.6)
C
      TOUT = DTOUT
      DO 70 JOUT = 1, 10
C
        CALL FCVODE(TOUT, T, U, ITASK, IER)
C
        CALL MAXNORM(NEQ, U, UNORM)
        WRITE(6,50) T, UNORM, IOUT(LNST)
 50     FORMAT(' At t = ', F6.2, '  max.norm(u) = ', E14.6,
     1         '  NST = ', I4)
C
        IF (IER .NE. 0) THEN
          WRITE(6,60) IER, IOUT(15)
 60       FORMAT(///' SUNDIALS_ERROR: FCVODE returned IER = ', I5, /,
     1           '                 Linear Solver returned IER = ', I5)
          CALL FCVFREE
          STOP
          ENDIF
C
          TOUT = TOUT + DTOUT
 70    CONTINUE
C
      WRITE(6,80) IOUT(LNST), IOUT(LNFE), IOUT(LNJE), IOUT(LNSETUP),
     1            IOUT(LNNI), IOUT(LNCF), IOUT(LNETF)
 80   FORMAT(//'Final statistics:'//
     1       ' No. steps = ', I4, '  No. f-s = ', I4,
     2       '  No. J-s = ', I4, '   No. LU-s = ', I4/
     3       ' No. nonlinear iterations = ', I4/
     4       ' No. nonlinear convergence failures = ', I4/
     5       ' No. error test failures = ', I4)
C
      CALL FCVFREE
C
      STOP
      END

C     ----------------------------------------------------------------

      SUBROUTINE INITBX(XMAX, YMAX, MX, MY, U0, IPAR, RPAR)
C Load IPAR and RPAR with problem constants and U0 with initial values
      IMPLICIT NONE
C
C The following declaration specification should match C type long int.
      INTEGER*8 IPAR(*)
      INTEGER*4 MX, MY
      DOUBLE PRECISION XMAX, YMAX, U0(MY,MX), RPAR(*)
C
      INTEGER*4 I, J
      DOUBLE PRECISION DX, DY, X, Y, HDCOEF, HACOEF, VDCOEF
C
C Problem constants
      DX = XMAX / (MX + 1)
      DY = YMAX / (MY + 1)
      HDCOEF = 1.0D0 / (DX * DX)
      HACOEF = 0.5D0 / (2.0D0 * DX)
      VDCOEF = 1.0D0 / (DY * DY)
C Load constants in IPAR and RPAR
      IPAR(1) = MX
      IPAR(2) = MY
      RPAR(1) = DX
      RPAR(2) = DY
      RPAR(3) = HDCOEF
      RPAR(4) = HACOEF
      RPAR(5) = VDCOEF
C
C Loop over grid and load initial values.
      DO 20 I = 1, MX
        X = I * DX
        DO 10 J = 1, MY
          Y = J * DY
          U0(J,I) = X * (XMAX - X) * Y * (YMAX - Y) *
     *              EXP(5.0D0 * X * Y)
 10       CONTINUE
 20     CONTINUE
C
      RETURN
      END

      SUBROUTINE MAXNORM(N, U, UNORM)
C Compute max-norm of array U
      IMPLICIT NONE
C
      INTEGER*8 I, N
      DOUBLE PRECISION U(*), UNORM, TEMP
C
      TEMP = 0.0D0
      DO 10 I = 1, N
         TEMP = MAX(ABS(U(I)), TEMP)
 10   CONTINUE
      UNORM = TEMP
      RETURN
      END

C     ----------------------------------------------------------------

      SUBROUTINE FCVFUN(T, U, UDOT, IPAR, RPAR, IER)
C Right-hand side routine
      IMPLICIT NONE
C
      DOUBLE PRECISION T, U(*), UDOT(*), RPAR(*)
C The following declaration specification should match C type long int.
      INTEGER*8 IPAR(*)
      INTEGER*4 IER
C
      INTEGER*4 I, MX, IOFF, MY, J, IJ
      DOUBLE PRECISION UIJ, UDN, UUP, ULT, URT, HDIFF, HADV, VDIFF
      DOUBLE PRECISION DX, DY, HDCOEF, HACOEF, VDCOEF
C
C Extract constants from IPAR and RPAR
      MX     = IPAR(1)
      MY     = IPAR(2)
      DX     = RPAR(1)
      DY     = RPAR(2)
      HDCOEF = RPAR(3)
      HACOEF = RPAR(4)
      VDCOEF = RPAR(5)
C
C Loop over all grid points.
      DO 20 I = 1, MX
        IOFF = (I - 1) * MY
        DO 10 J = 1, MY
C
C Extract u at x_i, y_j and four neighboring points.
          IJ = J + IOFF
          UIJ = U(IJ)
          UDN = 0.0D0
          IF (J .NE. 1)  UDN = U(IJ - 1)
          UUP = 0.0D0
          IF (J .NE. MY) UUP = U(IJ + 1)
          ULT = 0.0D0
          IF (I .NE. 1)  ULT = U(IJ - MY)
          URT = 0.0D0
          IF (I .NE. MX) URT = U(IJ + MY)
C
C Set diffusion and advection terms and load into UDOT.
          HDIFF = HDCOEF * (ULT - 2.0D0 * UIJ + URT)
          HADV = HACOEF * (URT - ULT)
          VDIFF = VDCOEF * (UUP - 2.0D0 * UIJ + UDN)
          UDOT(IJ) = HDIFF + HADV + VDIFF
 10       CONTINUE
 20     CONTINUE
C
        IER = 0
C
      RETURN
      END

C     ----------------------------------------------------------------

      SUBROUTINE FCVBJAC(N, MU, ML, MDIM, T, U, FU,
     1                   BJAC, H, IPAR, RPAR, V1, V2, V3, IER)
C Load banded Jacobian
      IMPLICIT NONE
C
C The following declaration specification should match C type long int.
      INTEGER*8 N, MU, ML, MDIM, IPAR(*)
      INTEGER*4 IER
      DOUBLE PRECISION T, U(*), FU(*), BJAC(MDIM,*), H, RPAR(*)
      DOUBLE PRECISION V1(*), V2(*), V3(*)
C
      INTEGER*4 MBAND, MX, MY
      INTEGER*4 I, J, K, IOFF, MU1, MU2
      DOUBLE PRECISION DX, DY, HDCOEF, HACOEF, VDCOEF
C
C Extract constants from IPAR and RPAR
      MX     = IPAR(1)
      MY     = IPAR(2)
      DX     = RPAR(1)
      DY     = RPAR(2)
      HDCOEF = RPAR(3)
      HACOEF = RPAR(4)
      VDCOEF = RPAR(5)
C
      MU1 = MU + 1
      MU2 = MU + 2
      MBAND = MU + 1 + ML
C
C Loop over all grid points.
      DO 20 I = 1, MX
        IOFF = (I - 1) * MY
        DO 10 J = 1, MY
          K = J + IOFF
C
C Set Jacobian elements in column k of Jb.
          BJAC(MU1,K) = -2.0D0 * (VDCOEF + HDCOEF)
          IF (I .NE. 1)  BJAC(1,K) = HDCOEF + HACOEF
          IF (I .NE. MX) BJAC(MBAND,K) = HDCOEF - HACOEF
          IF (J .NE. 1)  BJAC(MU,K) = VDCOEF
          IF (J .NE. MY) BJAC(MU2,K) = VDCOEF
C
 10       CONTINUE
 20     CONTINUE
C
        IER = 0
C     
      RETURN
      END
