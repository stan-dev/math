C     ----------------------------------------------------------------
C     Programmer(s): Daniel R. Reynolds @ SMU
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
C     ----------------------------------------------------------------
C     FARKODE Example Problem: 2D kinetics-transport, 
C     precond. Krylov solver. 
C     
C     An ODE system is generated from the following 2-species diurnal
C     kinetics advection-diffusion PDE system in 2 space dimensions:
C     
C     dc(i)/dt = Kh*(d/dx)**2 c(i) + V*dc(i)/dx + (d/dy)(Kv(y)*dc(i)/dy)
C                           + Ri(c1,c2,t)      for i = 1,2,   where
C     R1(c1,c2,t) = -q1*c1*c3 - q2*c1*c2 + 2*q3(t)*c3 + q4(t)*c2 ,
C     R2(c1,c2,t) =  q1*c1*c3 - q2*c1*c2 - q4(t)*c2 ,
C     Kv(y) = Kv0*exp(y/5) ,
C     Kh, V, Kv0, q1, q2, and c3 are constants, and q3(t) and q4(t)
C     vary diurnally.
C
C     The problem is posed on the square
C     0 .le. x .le. 20,    30 .le. y .le. 50   (all in km),
C     with homogeneous Neumann boundary conditions, and for time t in
C     0 .le. t .le. 86400 sec (1 day).
C
C     The PDE system is treated by central differences on a uniform
C     10 x 10 mesh, with simple polynomial initial profiles.
C     The problem is solved with ARKODE, with the DIRK/GMRES method and
C     using the FARKBP banded preconditioner.
C
C     Note that this problem should only work with SUNDIALS configured
C     to use 'realtype' as 'double' and 'sunindextype' as '64bit'
C     
C     The second and third dimensions of U here must match the values of
C     MX and MY, for consistency with the output statements below.
C     ----------------------------------------------------------------
C
      IMPLICIT NONE
      INTEGER*4 MX, MY
      PARAMETER (MX=10, MY=10)
      INTEGER*8 NEQ
      PARAMETER (NEQ=2*MX*MY)
C
      INTEGER*4 LNST, LNST_ATT, LNFE, LNFI, LNSETUP, LNNI, LNCF, LNPE
      INTEGER*4 LNLI, LNPS, LNCFL, LH, LNETF, LLENRW, LLENIW, LLENRWLS
      INTEGER*4 LLENIWLS, METH, IATOL, ITASK, IER, MAXL, JPRETYPE
      INTEGER*4 IGSTYPE, JOUT
      INTEGER*8 IOUT(35), IPAR(4), NST, NST_ATT, NFE, NFI, NPSET, NPE
      INTEGER*8 NPS, NNI, NLI, NCFN, NCFL, NETF, MU, ML, LENRW, LENIW
      INTEGER*8 LENRWLS, LENIWLS, LENRWBP, LENIWBP, NFEBP, MXSTEPS
      DOUBLE PRECISION ATOL, AVDIM, DELT, FLOOR, RTOL, T, TOUT, TWOHR
      DOUBLE PRECISION ROUT(6), U(2,MX,MY), RPAR(12)
C
      DATA TWOHR/7200.0D0/, RTOL/1.0D-5/, FLOOR/100.0D0/,
     1     JPRETYPE/1/, IGSTYPE/1/, MAXL/0/, DELT/0.0D0/, MXSTEPS/10000/
      DATA LLENRW/1/, LLENIW/2/, LNST/3/, LNST_ATT/6/, LNFE/7/, LNFI/8/, 
     1     LNETF/10/, LNCF/12/, LNNI/11/, LNSETUP/9/, LLENRWLS/14/, 
     1     LLENIWLS/15/, LNPE/21/, LNLI/23/, LNPS/22/, LNCFL/24/
      DATA LH/2/
C
C Load IPAR, RPAR, and initial values
      CALL INITKX(MX, MY, U, IPAR, RPAR)
C
C     Set other input arguments.
      T = 0.0D0
      METH = 0
      IATOL = 1
      ATOL = RTOL * FLOOR
      ITASK = 1
C
      WRITE(6,10) NEQ
 10   FORMAT('Krylov example problem:'//
     1       ' Kinetics-transport, NEQ = ', I4/)
C     
C     Initialize vector specification
      CALL FNVINITS(4, NEQ, IER)
      IF (IER .NE. 0) THEN
         WRITE(6,20) IER
 20      FORMAT(///' SUNDIALS_ERROR: FNVINITS returned IER = ', I5)
         STOP
      ENDIF
C
C     initialize SPGMR linear solver module
      call FSUNSPGMRINIT(4, JPRETYPE, MAXL, IER)
      IF (IER .NE. 0) THEN
        WRITE(6,25) IER
 25     FORMAT(///' SUNDIALS_ERROR: FSUNSPGMRINIT IER = ', I5)
        STOP
      ENDIF
      call FSUNSPGMRSETGSTYPE(4, IGSTYPE, IER)
      IF (IER .NE. 0) THEN
        WRITE(6,27) IER
 27     FORMAT(///' SUNDIALS_ERROR: FSUNSPGMRSETGSTYPE IER = ', I5)
        STOP
      ENDIF
C     
C     Initialize ARKODE
      CALL FARKMALLOC(T, U, METH, IATOL, RTOL, ATOL,
     1               IOUT, ROUT, IPAR, RPAR, IER)
      IF (IER .NE. 0) THEN
         WRITE(6,30) IER
 30      FORMAT(///' SUNDIALS_ERROR: FARKMALLOC returned IER = ', I5)
         STOP
      ENDIF
C
      CALL FARKSETIIN('MAX_NSTEPS', MXSTEPS, IER)
      IF (IER .NE. 0) THEN
        WRITE(6,35) IER
 35     FORMAT(///' SUNDIALS_ERROR: FARKSETIIN returned IER = ', I5)
        STOP
      ENDIF
C
C     attach matrix and linear solver modules to ARKLs interface
      CALL FARKLSINIT(IER)
      IF (IER .NE. 0) THEN
        WRITE(6,40) IER
 40     FORMAT(///' SUNDIALS_ERROR: FARKLSINIT returned IER = ',I5)
        CALL FARKFREE
        STOP
      ENDIF
C     
C     Initialize band preconditioner
      MU = 2
      ML = 2
      CALL FARKBPINIT(NEQ, MU, ML, IER) 
      IF (IER .NE. 0) THEN
         WRITE(6,45) IER
 45      FORMAT(///' SUNDIALS_ERROR: FARKBPINIT returned IER = ', I5)
         CALL FARKFREE
         STOP
      ENDIF
C     
C     Loop over output points, call FARKODE, print sample solution values.
      TOUT = TWOHR
      DO 70 JOUT = 1, 12
C
         CALL FARKODE(TOUT, T, U, ITASK, IER)
C     
         WRITE(6,50) T, IOUT(LNST), IOUT(LNST_ATT), ROUT(LH)
 50      FORMAT(/' t = ', E14.6, 5X, 'no. steps = ', I5, 
     1        '  no. att. steps = ', I5, '   stepsize = ', E14.6)
         WRITE(6,55) U(1,1,1), U(1,5,5), U(1,10,10),
     1               U(2,1,1), U(2,5,5), U(2,10,10)
 55      FORMAT('  c1 (bot.left/middle/top rt.) = ', 3E14.6/
     1        '  c2 (bot.left/middle/top rt.) = ', 3E14.6)
C     
         IF (IER .NE. 0) THEN
            WRITE(6,60) IER, IOUT(16)
 60         FORMAT(///' SUNDIALS_ERROR: FARKODE returned IER = ', I5, /,
     1             '                 Linear Solver returned IER = ', I5)
            CALL FARKFREE
            STOP
         ENDIF
C     
         TOUT = TOUT + TWOHR
 70   CONTINUE
      
C     Print final statistics.
      NST = IOUT(LNST)
      NST_ATT = IOUT(LNST_ATT)
      NFE = IOUT(LNFE)
      NFI = IOUT(LNFI)
      NPSET = IOUT(LNSETUP)
      NPE = IOUT(LNPE)
      NPS = IOUT(LNPS)
      NNI = IOUT(LNNI)
      NLI = IOUT(LNLI)
      AVDIM = DBLE(NLI) / DBLE(NNI)
      NCFN = IOUT(LNCF)
      NCFL = IOUT(LNCFL)
      NETF = IOUT(LNETF)
      LENRW = IOUT(LLENRW)
      LENIW = IOUT(LLENIW)
      LENRWLS = IOUT(LLENRWLS)
      LENIWLS = IOUT(LLENIWLS)
      WRITE(6,80) NST, NST_ATT, NFE, NFI, NPSET, NPE, NPS, NNI, NLI, 
     1     AVDIM, NCFN, NCFL, NETF, LENRW, LENIW, LENRWLS, LENIWLS
 80   FORMAT(//'Final statistics:'//
     &   ' number of steps        = ', I6/
     &   ' number of step att.    = ', I6/
     &   ' number of fe evals.    = ', I6/
     &   ' number of fi evals.    = ', I6/
     &   ' number of prec. setups = ', I6/
     &   ' number of prec. evals. = ', I6/
     &   ' number of prec. solves = ', I6/
     &   ' number of nonl. iters. = ', I6/
     &   ' number of lin. iters.  = ', I6/
     &   ' average Krylov subspace dimension (NLI/NNI) = ', E14.6/
     &   ' number of conv. failures.. nonlinear =', I3,
     &   ' linear = ', I3/
     &   ' number of error test failures = ', I3/
     &   ' main solver real/int workspace sizes   = ',2I5/
     &   ' linear solver real/int workspace sizes = ',2I5)
      CALL FARKBPOPT(LENRWBP, LENIWBP, NFEBP)
      WRITE(6,82) LENRWBP, LENIWBP, NFEBP
 82   FORMAT('In ARKBANDPRE:'/
     &        ' real/int workspace sizes = ', 2I5/
     &        ' number of f evaluations  = ', I5)
C     
      CALL FARKFREE
C     
      STOP
      END

C     ----------------------------------------------------------------

      SUBROUTINE INITKX(MX, MY, U0, IPAR, RPAR)
C Routine to set problem constants and initial values
C
      IMPLICIT NONE
C
      INTEGER*4 MX, MY
      INTEGER*8 IPAR(*), NEQ
      DOUBLE PRECISION RPAR(*)
C
      INTEGER*8 MM, JY, JX
      DOUBLE PRECISION U0
      DIMENSION U0(2,MX,MY)
      DOUBLE PRECISION Q1, Q2, Q3, Q4, A3, A4, OM, C3, DY, HDCO
      DOUBLE PRECISION VDCO, HACO, X, Y
      DOUBLE PRECISION CX, CY, DKH, DKV0, DX, HALFDA, PI, VEL
C
      DATA DKH/4.0D-6/, VEL/0.001D0/, DKV0/1.0D-8/, HALFDA/4.32D4/,
     1     PI/3.1415926535898D0/
C
C Problem constants
      MM = MX * MY
      NEQ = 2 * MM
      Q1 = 1.63D-16
      Q2 = 4.66D-16
      Q3 = 0.D0
      Q4 = 0.D0
      A3 = 22.62D0
      A4 = 7.601D0
      OM = PI / HALFDA
      C3 = 3.7D16
      DX = 20.0D0 / (MX - 1.0D0)
      DY = 20.0D0 / (MY - 1.0D0)
      HDCO = DKH / DX**2
      HACO = VEL / (2.0D0 * DX)
      VDCO = (1.0D0 / DY**2) * DKV0
C Load constants in IPAR and RPAR
      IPAR(1) = MX
      IPAR(2) = MY
      IPAR(3) = MM
      IPAR(4) = NEQ
C
      RPAR(1)  = Q1
      RPAR(2)  = Q2
      RPAR(3)  = Q3
      RPAR(4)  = Q4
      RPAR(5)  = A3
      RPAR(6)  = A4
      RPAR(7)  = OM
      RPAR(8)  = C3
      RPAR(9)  = DY
      RPAR(10) = HDCO
      RPAR(11) = VDCO
      RPAR(12) = HACO
C
C Set initial profiles.
      DO 20 JY = 1, MY
        Y = 30.0D0 + (JY - 1.0D0) * DY
        CY = (0.1D0 * (Y - 40.0D0))**2
        CY = 1.0D0 - CY + 0.5D0 * CY**2
        DO 10 JX = 1, MX
          X = (JX - 1.0D0) * DX
          CX = (0.1D0 * (X - 10.0D0))**2
          CX = 1.0D0 - CX + 0.5D0 * CX**2
          U0(1,JX,JY) = 1.0D6 * CX * CY
          U0(2,JX,JY) = 1.0D12 * CX * CY
 10       CONTINUE
 20     CONTINUE
C
      RETURN
      END
      
C     ----------------------------------------------------------------

      SUBROUTINE FARKIFUN(T, U, UDOT, IPAR, RPAR, IER)
C     Routine for right-hand side function fi
      IMPLICIT NONE
C
      INTEGER*4 IER
      INTEGER*8 IPAR(*)
      DOUBLE PRECISION T, U(2,*), UDOT(2,*), RPAR(*)
C
      INTEGER*4 ILEFT, IRIGHT
      INTEGER*8 MX, MY, MM, JY, JX, IBLOK0, IDN, IUP, IBLOK
      DOUBLE PRECISION Q1,Q2,Q3,Q4, A3, A4, OM, C3, DY, HDCO, VDCO, HACO
      DOUBLE PRECISION C1, C2, C1DN, C2DN, C1UP, C2UP, C1LT, C2LT
      DOUBLE PRECISION C1RT, C2RT, CYDN, CYUP, HORD1, HORD2, HORAD1
      DOUBLE PRECISION HORAD2, QQ1, QQ2, QQ3, QQ4, RKIN1, RKIN2, S
      DOUBLE PRECISION VERTD1, VERTD2, YDN, YUP
C
C     Extract constants from IPAR and RPAR
      MX = IPAR(1)
      MY = IPAR(2)
      MM = IPAR(3)
C
      Q1 = RPAR(1)
      Q2 = RPAR(2)
      Q3 = RPAR(3)
      Q4 = RPAR(4)
      A3 = RPAR(5)
      A4 = RPAR(6)
      OM = RPAR(7)
      C3 = RPAR(8)
      DY = RPAR(9)
      HDCO = RPAR(10)
      VDCO = RPAR(11)
      HACO = RPAR(12)
C     
C     Set diurnal rate coefficients.
      S = SIN(OM * T)
      IF (S .GT. 0.0D0) THEN
         Q3 = EXP(-A3 / S)
         Q4 = EXP(-A4 / S)
      ELSE
         Q3 = 0.0D0
         Q4 = 0.0D0
      ENDIF
C     
C     Loop over all grid points.
      DO 20 JY = 1, MY
         YDN = 30.0D0 + (JY - 1.5D0) * DY
         YUP = YDN + DY
         CYDN = VDCO * EXP(0.2D0 * YDN)
         CYUP = VDCO * EXP(0.2D0 * YUP)
         IBLOK0 = (JY - 1) * MX
         IDN = -MX
         IF (JY .EQ. 1) IDN = MX
         IUP = MX
         IF (JY .EQ. MY) IUP = -MX
         DO 10 JX = 1, MX
            IBLOK = IBLOK0 + JX
            C1 = U(1,IBLOK)
            C2 = U(2,IBLOK)
C     Set kinetic rate terms.
            QQ1 = Q1 * C1 * C3
            QQ2 = Q2 * C1 * C2
            QQ3 = Q3 * C3
            QQ4 = Q4 * C2
            RKIN1 = -QQ1 - QQ2 + 2.0D0 * QQ3 + QQ4
            RKIN2 = QQ1 - QQ2 - QQ4
C     Set vertical diffusion terms.
            C1DN = U(1,IBLOK + IDN)
            C2DN = U(2,IBLOK + IDN)
            C1UP = U(1,IBLOK + IUP)
            C2UP = U(2,IBLOK + IUP)
            VERTD1 = CYUP * (C1UP - C1) - CYDN * (C1 - C1DN)
            VERTD2 = CYUP * (C2UP - C2) - CYDN * (C2 - C2DN)
C     Set horizontal diffusion and advection terms.
            ILEFT = -1
            IF (JX .EQ. 1) ILEFT = 1
            IRIGHT = 1
            IF (JX .EQ. MX) IRIGHT = -1
            C1LT = U(1,IBLOK + ILEFT)
            C2LT = U(2,IBLOK + ILEFT)
            C1RT = U(1,IBLOK + IRIGHT)
            C2RT = U(2,IBLOK + IRIGHT)
            HORD1 = HDCO * (C1RT - 2.0D0 * C1 + C1LT)
            HORD2 = HDCO * (C2RT - 2.0D0 * C2 + C2LT)
            HORAD1 = HACO * (C1RT - C1LT)
            HORAD2 = HACO * (C2RT - C2LT)
C     Load all terms into UDOT.
            UDOT(1,IBLOK) = VERTD1 + HORD1 + HORAD1 + RKIN1
            UDOT(2,IBLOK) = VERTD2 + HORD2 + HORAD2 + RKIN2
 10      CONTINUE
 20   CONTINUE
C
      IER = 0
C
      RETURN
      END
      
C     ----------------------------------------------------------------

      SUBROUTINE FARKEFUN(T, U, UDOT, IPAR, RPAR, IER)
C     Routine for right-hand side function fe
      IMPLICIT NONE
C
      INTEGER*4 IER
      INTEGER*8 IPAR(*)
      DOUBLE PRECISION T, U(2,*), UDOT(2,*), RPAR(*)
C
      INTEGER*8 MX, MY, JY, JX, IBLOK0, IBLOK
C
C     Extract constants from IPAR
      MX = IPAR(1)
      MY = IPAR(2)
C     
C     Loop over all grid points, setting UDOT to zero
      DO 20 JY = 1, MY
         IBLOK0 = (JY - 1) * MX
         DO 10 JX = 1, MX
            IBLOK = IBLOK0 + JX
            UDOT(1,IBLOK) = 0.D0
            UDOT(2,IBLOK) = 0.D0
 10      CONTINUE
 20   CONTINUE
C
      IER = 0
C
      RETURN
      END
