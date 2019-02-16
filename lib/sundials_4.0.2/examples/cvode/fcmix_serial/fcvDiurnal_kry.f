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
C     FCVODE Example Problem: 2D kinetics-transport, precond. Krylov
C     solver. 
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
C     with homogeneous Neumann boundary conditions, and for time t
C     in 0 .le. t .le. 86400 sec (1 day).
C
C     The PDE system is treated by central differences on a uniform
C     10 x 10 mesh, with simple polynomial initial profiles.
C     The problem is solved with CVODE, with the BDF/GMRES method and
C     the block-diagonal part of the Jacobian as a left
C     preconditioner.
C     
C     Note: this program includes modified (64-bit integer) versions 
C     of the dense linear solver routines DGEFA and DGESL from LINPACK, 
C     and BLAS routines DCOPY and DSCAL.
C     
C     The second and third dimensions of U here must match the values
C     of MX and MY, for consistency with the output statements below.
C     --------------------------------------------------------------------
C
      IMPLICIT NONE
C
      INTEGER*4 MX, MY
      PARAMETER (MX=10, MY=10)
      INTEGER*4 LENIPAR, LENRPAR
      PARAMETER (LENIPAR=6+2*MX*MY, LENRPAR=12+8*MX*MY)
C
      INTEGER*4 METH,IATOL,ITASK,IER,LNCFL,LNPS
      INTEGER*4 LNST,LNFE,LNSETUP,LNNI,LNCF,LQ,LH,LNPE,LNLI,LNETF
      INTEGER*4 JOUT,JPRETYPE,IGSTYPE,MAXL
C The following declaration specification should match C type long int.
      INTEGER*8 NEQ, IOUT(25), IPAR(LENIPAR)
      INTEGER*4 NST,NFE,NPSET,NPE,NPS,NNI,NETF
      INTEGER*4 NLI,NCFN,NCFL
      DOUBLE PRECISION ATOL,AVDIM,T,TOUT,TWOHR,RTOL,FLOOR,DELT
      DOUBLE PRECISION U(2,MX,MY),ROUT(10),RPAR(LENRPAR)
C
      DATA TWOHR/7200.0D0/, RTOL/1.0D-5/, FLOOR/100.0D0/,
     &     JPRETYPE/1/, IGSTYPE/1/, MAXL/0/, DELT/0.0D0/
      DATA LNST/3/, LNFE/4/, LNETF/5/,  LNCF/6/, LNNI/7/, LNSETUP/8/, 
     &     LQ/9/, LNPE/20/, LNLI/22/, LNPS/21/, LNCFL/23/
      DATA LH/2/
C
C     Load problem constants into IPAR, RPAR, and set initial values
      CALL INITKX(MX, MY, U, IPAR, RPAR)
C
C     Set other input arguments.
      NEQ = 2*MX*MY
      T = 0.0D0
      METH = 2
      IATOL = 1
      ATOL = RTOL * FLOOR
      ITASK = 1
C
      WRITE(6,10) NEQ
 10   FORMAT('Krylov example problem:'//
     &       ' Kinetics-transport, NEQ = ', I4/)
C
C     Initialize vector specification
      CALL FNVINITS(1, NEQ, IER)
      IF (IER .NE. 0) THEN
        WRITE(6,20) IER
 20     FORMAT(///' SUNDIALS_ERROR: FNVINITS returned IER = ', I5)
        STOP
      ENDIF
C
C     Initialize SPGMR linear solver module
      call FSUNSPGMRINIT(1, JPRETYPE, MAXL, IER)
      IF (IER .NE. 0) THEN
        WRITE(6,25) IER
 25     FORMAT(///' SUNDIALS_ERROR: FSUNSPGMRINIT IER = ', I5)
        STOP
      ENDIF
      call FSUNSPGMRSETGSTYPE(1, IGSTYPE, IER)
      IF (IER .NE. 0) THEN
        WRITE(6,27) IER
 27     FORMAT(///' SUNDIALS_ERROR: FSUNSPGMRSETGSTYPE IER = ', I5)
        STOP
      ENDIF
C     
C     Initialize CVODE
      CALL FCVMALLOC(T, U, METH, IATOL, RTOL, ATOL,
     &     IOUT, ROUT, IPAR, RPAR, IER)
      IF (IER .NE. 0) THEN
        WRITE(6,30) IER
 30     FORMAT(///' SUNDIALS_ERROR: FCVMALLOC returned IER = ', I5)
        STOP
      ENDIF
C
C     attach linear solver module to CVLs interface
      CALL FCVLSINIT(IER)
      IF (IER .NE. 0) THEN
        WRITE(6,40) IER
 40     FORMAT(///' SUNDIALS_ERROR: FCVLSINIT returned IER = ',I5)
        CALL FCVFREE
        STOP
      ENDIF
C     
C     attach preconditioner to CVLs interface
      CALL FCVLSSETPREC(1, IER)
C
C Loop over output points, call FCVODE, print sample solution values.
      TOUT = TWOHR
      DO JOUT = 1, 12
C
         CALL FCVODE(TOUT, T, U, ITASK, IER)
C
         WRITE(6,50) T, IOUT(LNST), IOUT(LQ), ROUT(LH)
 50      FORMAT(/' t = ', E11.3, 3X, 'nst = ', I5,
     &           '  q = ', I2, '  h = ', E14.6)
         WRITE(6,55) U(1,1,1), U(1,5,5), U(1,10,10),
     &               U(2,1,1), U(2,5,5), U(2,10,10)
 55      FORMAT('  c1 (bot.left/middle/top rt.) = ', 3E14.6/
     &          '  c2 (bot.left/middle/top rt.) = ', 3E14.6)
C
         IF (IER .NE. 0) THEN
            WRITE(6,60) IER, IOUT(15)
 60         FORMAT(///' SUNDIALS_ERROR: FCVODE returned IER = ', I5, /,
     &             '                 Linear Solver returned IER = ', I5)
            CALL FCVFREE
            STOP
         ENDIF
C     
         TOUT = TOUT + TWOHR
C
      ENDDO

C     Print final statistics.
      NST = IOUT(LNST)
      NFE = IOUT(LNFE)
      NPSET = IOUT(LNSETUP)
      NPE = IOUT(LNPE)
      NPS = IOUT(LNPS)
      NNI = IOUT(LNNI)
      NLI = IOUT(LNLI)
      AVDIM = DBLE(NLI) / DBLE(NNI)
      NCFN = IOUT(LNCF)
      NCFL = IOUT(LNCFL)
      NETF = IOUT(LNETF)
      WRITE(6,80) NST, NFE, NPSET, NPE, NPS, NNI, NLI, AVDIM, NCFN,
     &     NCFL, NETF
  80  FORMAT(//'Final statistics:'//
     &     ' number of steps        = ', I5, 5X,
     &     ' number of f evals.     =', I5/
     &     ' number of prec. setups = ', I5/
     &     ' number of prec. evals. = ', I5, 5X,
     &     ' number of prec. solves = ', I5/
     &     ' number of nonl. iters. = ', I5, 5X,
     &     ' number of lin. iters.  = ', I5/
     &     ' average Krylov subspace dimension (NLI/NNI)  = ', E14.6/
     &     ' number of conv. failures.. nonlinear = ', I3,
     &     ' linear = ', I3/
     &     ' number of error test failures = ', I3)
C     
      CALL FCVFREE
C
      STOP
      END

C     ----------------------------------------------------------------

      SUBROUTINE INITKX(MX, MY, U0, IPAR, RPAR)
C     Routine to set problem constants and initial values
C
      IMPLICIT NONE
C
      INTEGER*4 MX, MY
C The following declaration specification should match C type long int.
      INTEGER*8 IPAR(*)
      DOUBLE PRECISION RPAR(*)
C
      INTEGER*4 MM, JY, JX, P_IPP, P_BD, P_P
      DOUBLE PRECISION U0
      DIMENSION U0(2,MX,MY)
      DOUBLE PRECISION Q1, Q2, Q3, Q4, A3, A4, OM, C3, DY, HDCO
      DOUBLE PRECISION VDCO, HACO, X, Y
      DOUBLE PRECISION CX, CY, DKH, DKV0, DX, HALFDA, PI, VEL
C
      DATA DKH/4.0D-6/, VEL/0.001D0/, DKV0/1.0D-8/, HALFDA/4.32D4/,
     1     PI/3.1415926535898D0/
C
C     Problem constants
      MM = MX * MY
      Q1 = 1.63D-16
      Q2 = 4.66D-16
      Q3 = 0.0D0
      Q4 = 0.0D0
      A3 = 22.62D0
      A4 = 7.601D0
      OM = PI / HALFDA
      C3 = 3.7D16
      DX = 20.0D0 / (MX - 1.0D0)
      DY = 20.0D0 / (MY - 1.0D0)
      HDCO = DKH / DX**2
      HACO = VEL / (2.0D0 * DX)
      VDCO = (1.0D0 / DY**2) * DKV0
C
C     Load constants in IPAR and RPAR
      IPAR(1) = MX
      IPAR(2) = MY
      IPAR(3) = MM
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
C     Pointers into IPAR and RPAR
      P_IPP = 7
      P_BD  = 13
      P_P   = P_BD + 4*MM
C     
      IPAR(4) = P_IPP
      IPAR(5) = P_BD
      IPAR(6) = P_P
C
C     Set initial profiles.
      DO JY = 1, MY
        Y = 30.0D0 + (JY - 1.0D0) * DY
        CY = (0.1D0 * (Y - 40.0D0))**2
        CY = 1.0D0 - CY + 0.5D0 * CY**2
        DO JX = 1, MX
          X = (JX - 1.0D0) * DX
          CX = (0.1D0 * (X - 10.0D0))**2
          CX = 1.0D0 - CX + 0.5D0 * CX**2
          U0(1,JX,JY) = 1.0D6 * CX * CY
          U0(2,JX,JY) = 1.0D12 * CX * CY
         ENDDO
      ENDDO
C
      RETURN
      END

C     ----------------------------------------------------------------

      SUBROUTINE FCVFUN(T, U, UDOT, IPAR, RPAR, IER)
C     Routine for right-hand side function f
C
      IMPLICIT NONE
C
      DOUBLE PRECISION T, U(2,*), UDOT(2,*), RPAR(*)
C The following declaration specification should match C type long int.
      INTEGER*8 IPAR(*)
      INTEGER*4 IER
C
      INTEGER*4 ILEFT, IRIGHT
      INTEGER*4 JX, JY, MX, MY, MM, IBLOK0, IBLOK, IDN, IUP
      DOUBLE PRECISION Q1, Q2, Q3, Q4, A3, A4, OM, C3, DY, HDCO
      DOUBLE PRECISION VDCO, HACO
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
      RPAR(3) = Q3
      RPAR(4) = Q4
C
C     Loop over all grid points.
      DO JY = 1, MY
         YDN = 30.0D0 + (JY - 1.5D0) * DY
         YUP = YDN + DY
         CYDN = VDCO * EXP(0.2D0 * YDN)
         CYUP = VDCO * EXP(0.2D0 * YUP)
         IBLOK0 = (JY - 1) * MX
         IDN = -MX
         IF (JY .EQ. 1) IDN = MX
         IUP = MX
         IF (JY .EQ. MY) IUP = -MX
         DO JX = 1, MX
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
         ENDDO
      ENDDO
C
      IER = 0
C
      RETURN
      END

C     ----------------------------------------------------------------

      SUBROUTINE FCVPSET(T, U, FU, JOK, JCUR, GAMMA, H,
     &                   IPAR, RPAR, IER)
C     Routine to set and preprocess block-diagonal preconditioner.
C     Note: The dimensions in /BDJ/ below assume at most 100 mesh points.
C
      IMPLICIT NONE
C
      INTEGER*4 IER, JOK, JCUR
      DOUBLE PRECISION T, U(2,*), FU(*), GAMMA, H
C The following declaration specification should match C type long int.
      INTEGER*8 IPAR(*)
      DOUBLE PRECISION RPAR(*)
C
      INTEGER*4 MX, MY, MM, P_IPP, P_BD, P_P
      DOUBLE PRECISION Q1, Q2, Q3, Q4, C3, DY, HDCO, VDCO
C
      IER = 0
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
      C3 = RPAR(8)
      DY = RPAR(9)
      HDCO = RPAR(10)
      VDCO = RPAR(11)
C
C     Extract pointers into IPAR and RPAR
      P_IPP = IPAR(4)
      P_BD  = IPAR(5)
      P_P   = IPAR(6)
C
C     If needed, recompute BD
C
      IF (JOK .EQ. 1) THEN
C     JOK = 1. Reuse saved BD
        JCUR = 0
      ELSE    
C     JOK = 0. Compute diagonal Jacobian blocks.
C     (using q4 value computed on last FCVFUN call).
         CALL PREC_JAC(MX, MY, MM, U, RPAR(P_BD),
     &        Q1, Q2, Q3, Q4, C3, DY, HDCO, VDCO)
         JCUR = 1
      ENDIF
C
C     Copy BD to P
      CALL DCOPY(4*MM, RPAR(P_BD), 1, RPAR(P_P), 1)
C
C     Scale P by -GAMMA
      CALL DSCAL(4*MM, -GAMMA, RPAR(P_P), 1)
C
C     Perform LU decomposition
      CALL PREC_LU(MM, RPAR(P_P), IPAR(P_IPP), IER)
C
      RETURN
      END

C     ----------------------------------------------------------------
      
      SUBROUTINE FCVPSOL(T, U, FU, R, Z, GAMMA, DELTA, LR,
     &                   IPAR, RPAR, IER)
C     Routine to solve preconditioner linear system.
C
      IMPLICIT NONE
C
      INTEGER*4 IER, LR
C The following declaration specification should match C type long int.
      INTEGER*8 IPAR(*)
      DOUBLE PRECISION T, U(*), FU(*), R(*), Z(2,*)
      DOUBLE PRECISION GAMMA, DELTA, RPAR(*)
C
      INTEGER*4 MM, P_IPP, P_P
C
      IER = 0
C
C     Extract constants from IPAR and RPAR
      MM  = IPAR(3)
C
C     Extract pointers into IPAR and RPAR
      P_IPP = IPAR(4)
      P_P   = IPAR(6)
C
C     Copy RHS into Z
      CALL DCOPY(2*MM, R, 1, Z, 1)
C
C     Solve the block-diagonal system Px = r using LU factors stored in P
C     and pivot data in IPP, and return the solution in Z.
      CALL PREC_SOL(MM, RPAR(P_P), IPAR(P_IPP), Z)

      RETURN
      END

C     ----------------------------------------------------------------

      SUBROUTINE PREC_JAC(MX, MY, MM, U, BD,
     &     Q1, Q2, Q3, Q4, C3, DY, HDCO, VDCO)
C     Routine to compute diagonal Jacobian blocks
C
      IMPLICIT NONE
C
      INTEGER*4 MX, MY, MM
      DOUBLE PRECISION U(2,*), BD(2,2,MM)
      DOUBLE PRECISION Q1, Q2, Q3, Q4, C3, DY, HDCO, VDCO
C
      INTEGER*4 JY, JX, IBLOK, IBLOK0
      DOUBLE PRECISION C1, C2, CYDN, CYUP, DIAG, YDN, YUP
C
      DO JY = 1, MY
         YDN = 30.0D0 + (JY - 1.5D0) * DY
         YUP = YDN + DY
         CYDN = VDCO * EXP(0.2D0 * YDN)
         CYUP = VDCO * EXP(0.2D0 * YUP)
         DIAG = -(CYDN + CYUP + 2.0D0 * HDCO)
         IBLOK0 = (JY - 1) * MX
         DO JX = 1, MX
            IBLOK = IBLOK0 + JX
            C1 = U(1,IBLOK)
            C2 = U(2,IBLOK)
            BD(1,1,IBLOK) = (-Q1 * C3 - Q2 * C2) + DIAG
            BD(1,2,IBLOK) = -Q2 * C1 + Q4
            BD(2,1,IBLOK) =  Q1 * C3 - Q2 * C2
            BD(2,2,IBLOK) = (-Q2 * C1 - Q4) + DIAG
         ENDDO
      ENDDO

      RETURN
      END

C     ----------------------------------------------------------------

      SUBROUTINE PREC_LU(MM, P, IPP, IER)
C     Routine to perform LU decomposition on (P+I)
C
      IMPLICIT NONE
C
      INTEGER*4 IER
      INTEGER*4 MM
      INTEGER*8 IPP(2,MM)
      DOUBLE PRECISION P(2,2,MM)
C
      INTEGER*4 I
C
C     Add identity matrix and do LU decompositions on blocks, in place.
      DO I = 1, MM
         P(1,1,I) = P(1,1,I) + 1.0D0
         P(2,2,I) = P(2,2,I) + 1.0D0
         CALL DGEFA(P(1,1,I), 2, 2, IPP(1,I), IER)
         IF (IER .NE. 0) RETURN
      ENDDO
C     
      RETURN
      END

C     ----------------------------------------------------------------

      SUBROUTINE PREC_SOL(MM, P, IPP, Z)
C     Routine for backsolve
C
      IMPLICIT NONE
C
      INTEGER*4 MM
      INTEGER*8 IPP(2,MM)
      DOUBLE PRECISION P(2,2,MM), Z(2,MM)
C      
      INTEGER*4 I
C
      DO I = 1, MM
         CALL DGESL(P(1,1,I), 2, 2, IPP(1,I), Z(1,I), 0)
      ENDDO

      RETURN
      END

C     ----------------------------------------------------------------

      subroutine dgefa(a, lda, n, ipvt, info)
c
      implicit none
c
      integer info, idamax, j, k, kp1, l, nm1, n
      integer*4 lda
      integer*8 ipvt(1)
      double precision a(lda,1), t
c
c     dgefa factors a double precision matrix by gaussian elimination.
c
c     dgefa is usually called by dgeco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c     (time for dgeco) = (1 + 9/n)*(time for dgefa) .
c
c     on entry
c
c        a       double precision(lda, n)
c                the matrix to be factored.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that dgesl or dgedi will divide by zero
c                     if called.  employ  rcond  in dgeco for a reliable
c                     indication of singularity.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,dscal,idamax
c
c     internal variables
c
c     gaussian elimination with partial pivoting
c
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1
c
c        find l = pivot index
c
         l = idamax(n - k + 1, a(k,k), 1) + k - 1
         ipvt(k) = l
c
c        zero pivot implies this column already triangularized
c
         if (a(l,k) .eq. 0.0d0) go to 40
c
c           interchange if necessary
c
            if (l .eq. k) go to 10
               t = a(l,k)
               a(l,k) = a(k,k)
               a(k,k) = t
   10       continue
c
c           compute multipliers
c
            t = -1.0d0 / a(k,k)
            call dscal(n - k, t, a(k + 1,k), 1)
c
c           row elimination with column indexing
c
            do 30 j = kp1, n
               t = a(l,j)
               if (l .eq. k) go to 20
                  a(l,j) = a(k,j)
                  a(k,j) = t
   20          continue
               call daxpy(n - k, t, a(k + 1,k), 1, a(k + 1,j), 1)
   30       continue
         go to 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (a(n,n) .eq. 0.0d0) info = n
      return
      end

C     ----------------------------------------------------------------

      subroutine dgesl(a, lda, n, ipvt, b, job)
c
      implicit none
c
      integer lda, n, job, k, kb, l, nm1
      integer*8 ipvt(1)
      double precision a(lda,1), b(1), ddot, t
c
c     dgesl solves the double precision system
c     a * x = b  or  trans(a) * x = b
c     using the factors computed by dgeco or dgefa.
c
c     on entry
c
c        a       double precision(lda, n)
c                the output from dgeco or dgefa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        ipvt    integer(n)
c                the pivot vector from dgeco or dgefa.
c
c        b       double precision(n)
c                the right hand side vector.
c
c        job     integer
c                = 0         to solve  a*x = b ,
c                = nonzero   to solve  trans(a)*x = b  where
c                            trans(a)  is the transpose.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero will occur if the input factor contains a
c        zero on the diagonal.  technically this indicates singularity
c        but it is often caused by improper arguments or improper
c        setting of lda .  it will not occur if the subroutines are
c        called correctly and if dgeco has set rcond .gt. 0.0
c        or dgefa has set info .eq. 0 .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call dgeco(a,lda,n,ipvt,rcond,z)
c           if (rcond is too small) go to ...
c           do 10 j = 1, p
c              call dgesl(a,lda,n,ipvt,c(1,j),0)
c        10 continue
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,ddot
c
c     internal variables
c
      nm1 = n - 1
      if (job .ne. 0) go to 50
c
c        job = 0 , solve  a * x = b
c        first solve  l*y = b
c
         if (nm1 .lt. 1) go to 30
         do 20 k = 1, nm1
            l = ipvt(k)
            t = b(l)
            if (l .eq. k) go to 10
               b(l) = b(k)
               b(k) = t
   10       continue
            call daxpy(n - k, t, a(k + 1,k), 1, b(k + 1), 1)
   20    continue
   30    continue
c
c        now solve  u*x = y
c
         do 40 kb = 1, n
            k = n + 1 - kb
            b(k) = b(k) / a(k,k)
            t = -b(k)
            call daxpy(k - 1, t, a(1,k), 1, b(1), 1)
   40    continue
      go to 100
   50 continue
c
c        job = nonzero, solve  trans(a) * x = b
c        first solve  trans(u)*y = b
c
         do 60 k = 1, n
            t = ddot(k - 1, a(1,k), 1, b(1), 1)
            b(k) = (b(k) - t) / a(k,k)
   60    continue
c
c        now solve trans(l)*x = y
c
         if (nm1 .lt. 1) go to 90
         do 80 kb = 1, nm1
            k = n - kb
            b(k) = b(k) + ddot(n - k, a(k + 1,k), 1, b(k + 1), 1)
            l = ipvt(k)
            if (l .eq. k) go to 70
               t = b(l)
               b(l) = b(k)
               b(k) = t
   70       continue
   80    continue
   90    continue
  100 continue
      return
      end

C     ----------------------------------------------------------------

      subroutine daxpy(n, da, dx, incx, dy, incy)
c
c     constant times a vector plus a vector.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      implicit none
c
      integer i, incx, incy, ix, iy, m, mp1
      integer*4 n
      double precision dx(1), dy(1), da
c
      if (n .le. 0) return
      if (da .eq. 0.0d0) return
      if (incx .eq. 1 .and. incy .eq. 1) go to 20
c
c        code for unequal increments or equal increments
c        not equal to 1
c
      ix = 1
      iy = 1
      if (incx .lt. 0) ix = (-n + 1) * incx + 1
      if (incy .lt. 0) iy = (-n + 1) * incy + 1
      do 10 i = 1, n
        dy(iy) = dy(iy) + da * dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n, 4)
      if ( m .eq. 0 ) go to 40
      do 30 i = 1, m
        dy(i) = dy(i) + da * dx(i)
   30 continue
      if ( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1, n, 4
        dy(i) = dy(i) + da * dx(i)
        dy(i + 1) = dy(i + 1) + da * dx(i + 1)
        dy(i + 2) = dy(i + 2) + da * dx(i + 2)
        dy(i + 3) = dy(i + 3) + da * dx(i + 3)
   50 continue
      return
      end
C     ----------------------------------------------------------------

      subroutine dscal(n, da, dx, incx)
c
c     scales a vector by a constant.
c     uses unrolled loops for increment equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      implicit none
c
      integer i, incx, m, mp1, nincx
      integer*4 n
      double precision da, dx(1)
c
      if (n.le.0) return
      if (incx .eq. 1) go to 20
c
c        code for increment not equal to 1
c
      nincx = n * incx
      do 10 i = 1, nincx, incx
        dx(i) = da * dx(i)
   10 continue
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n, 5)
      if ( m .eq. 0 ) go to 40
      do 30 i = 1, m
        dx(i) = da * dx(i)
   30 continue
      if ( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1, n, 5
        dx(i) = da * dx(i)
        dx(i + 1) = da * dx(i + 1)
        dx(i + 2) = da * dx(i + 2)
        dx(i + 3) = da * dx(i + 3)
        dx(i + 4) = da * dx(i + 4)
   50 continue
      return
      end

C     ----------------------------------------------------------------

      double precision function ddot(n, dx, incx, dy, incy)
c
c     forms the dot product of two vectors.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      implicit none
c
      integer i, incx, incy, ix, iy, m, mp1
      integer*4 n
      double precision dx(1), dy(1), dtemp
c
      ddot = 0.0d0
      dtemp = 0.0d0
      if (n .le. 0) return
      if (incx .eq. 1 .and. incy .eq. 1) go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if (incx .lt. 0) ix = (-n + 1) * incx + 1
      if (incy .lt. 0) iy = (-n + 1) * incy + 1
      do 10 i = 1, n
        dtemp = dtemp + dx(ix) * dy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      ddot = dtemp
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n, 5)
      if ( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dtemp + dx(i) * dy(i)
   30 continue
      if ( n .lt. 5 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1, n, 5
        dtemp = dtemp + dx(i) * dy(i) + dx(i + 1) * dy(i + 1) +
     *          dx(i + 2) * dy(i + 2) + dx(i + 3) * dy(i + 3) +
     *          dx(i + 4) * dy(i + 4)
   50 continue
   60 ddot = dtemp
      return
      end

C     ----------------------------------------------------------------

      integer function idamax(n, dx, incx)
c
c     finds the index of element having max. absolute value.
c     jack dongarra, linpack, 3/11/78.
c
      implicit none
c
      integer i, incx, ix
      integer*4 n
      double precision dx(1), dmax
c
      idamax = 0
      if (n .lt. 1) return
      idamax = 1
      if (n .eq. 1) return
      if (incx .eq. 1) go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      dmax = abs(dx(1))
      ix = ix + incx
      do 10 i = 2, n
         if (abs(dx(ix)) .le. dmax) go to 5
         idamax = i
         dmax = abs(dx(ix))
    5    ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 dmax = abs(dx(1))
      do 30 i = 2, n
         if (abs(dx(i)) .le. dmax) go to 30
         idamax = i
         dmax = abs(dx(i))
   30 continue
      return
      end

C     ----------------------------------------------------------------

      subroutine  dcopy(n, dx, incx, dy, incy)
c
c     copies a vector, x, to a vector, y.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      implicit none
c
      integer i, incx, incy, ix, iy, m, mp1
      integer*4 n
      double precision dx(1), dy(1)
c
      if (n .le. 0) return
      if (incx .eq. 1 .and. incy .eq. 1) go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if (incx .lt. 0) ix = (-n + 1) * incx + 1
      if (incy .lt. 0) iy = (-n + 1) * incy + 1
      do 10 i = 1, n
        dy(iy) = dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n, 7)
      if ( m .eq. 0 ) go to 40
      do 30 i = 1, m
        dy(i) = dx(i)
   30 continue
      if ( n .lt. 7 ) return
   40 mp1 = m + 1
      do 50 i = mp1, n, 7
        dy(i) = dx(i)
        dy(i + 1) = dx(i + 1)
        dy(i + 2) = dx(i + 2)
        dy(i + 3) = dx(i + 3)
        dy(i + 4) = dx(i + 4)
        dy(i + 5) = dx(i + 5)
        dy(i + 6) = dx(i + 6)
   50 continue
      return
      end
