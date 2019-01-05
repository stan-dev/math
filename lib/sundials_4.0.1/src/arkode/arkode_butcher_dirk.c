/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
 * LLNS/SMU Copyright Start
 * Copyright (c) 2015, Southern Methodist University and
 * Lawrence Livermore National Security
 *
 * This work was performed under the auspices of the U.S. Department
 * of Energy by Southern Methodist University and Lawrence Livermore
 * National Laboratory under Contract DE-AC52-07NA27344.
 * Produced at Southern Methodist University and the Lawrence
 * Livermore National Laboratory.
 *
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS/SMU Copyright End
 *---------------------------------------------------------------
 * This is the implementation file for built-in DIRK Butcher
 * tables.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "arkode_impl.h"
#include <arkode/arkode_butcher_dirk.h>
#include <sundials/sundials_math.h>

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define RSYM ".32Lg"
#else
#define RSYM ".16g"
#endif


/*---------------------------------------------------------------
  Returns Butcher table structure for pre-set DIRK methods.

  Input:  imeth -- integer key for the desired method (see below)

  Allowed 'method' names and properties (those in an ARK pair are
  marked with a *).  All method names are of the form
  <name>_s_p_q.  The method 'type' is one of
    SDIRK -- singly-diagonally implicit Runge Kutta
    ESDIRK -- explicit [1st stage] singly-diagonally implicit Runge Kutta
  The 'A-stable' and 'L-stable' columns are based on numerical estimates
  of each property.  The 'QP' column denotes whether the coefficients
  of the method are known precisely enough for use in 'long double'
  (128-bit) calculations.

     imeth                       type  A-stable  L-stable  QP
    ----------------------------------------------------------
     SDIRK_2_1_2                SDIRK     Y         N       Y
     BILLINGTON_3_3_2           SDIRK     N         N       N
     TRBDF2_3_3_2              ESDIRK     N         N       Y
     KVAERNO_4_2_3             ESDIRK     Y         Y       N
     ARK324L2SA_DIRK_4_2_3*    ESDIRK     Y         Y       N
     CASH_5_2_4                 SDIRK     Y         Y       N
     CASH_5_3_4                 SDIRK     Y         Y       N
     SDIRK_5_3_4                SDIRK     Y         Y       Y
     KVAERNO_5_3_4             ESDIRK     Y         N       N
     ARK436L2SA_DIRK_6_3_4*    ESDIRK     Y         Y       N
     KVAERNO_7_4_5             ESDIRK     Y         Y       N
     ARK548L2SA_DIRK_8_4_5*    ESDIRK     Y         Y       N
    ----------------------------------------------------------

  ---------------------------------------------------------------*/
ARKodeButcherTable ARKodeButcherTable_LoadDIRK(int imethod)
{

  ARKodeButcherTable B;
  B = NULL;

  /* fill in coefficients based on method name */
  switch(imethod) {

  case(SDIRK_2_1_2):   /* SDIRK-2-1 (A,B stable) */
    B = ARKodeButcherTable_Alloc(2, SUNTRUE);
    B->q = 2;
    B->p = 1;

    B->A[0][0] = RCONST(1.0);
    B->A[1][0] = RCONST(-1.0);
    B->A[1][1] = RCONST(1.0);

    B->b[0] = RCONST(0.5);
    B->b[1] = RCONST(0.5);

    B->d[0] = RCONST(1.0);

    B->c[0] = RCONST(1.0);
    B->c[1] = RCONST(0.0);
    break;

  case(BILLINGTON_3_3_2):    /* Billington-SDIRK */
    B = ARKodeButcherTable_Alloc(3, SUNTRUE);
    B->q = 2;
    B->p = 3;

    B->A[0][0] = RCONST(0.292893218813);
    B->A[1][0] = RCONST(0.798989873223);
    B->A[1][1] = RCONST(0.292893218813);
    B->A[2][0] = RCONST(0.740789228841);
    B->A[2][1] = RCONST(0.259210771159);
    B->A[2][2] = RCONST(0.292893218813);

    B->d[0] = RCONST(0.691665115992);
    B->d[1] = RCONST(0.503597029883);
    B->d[2] = RCONST(-0.195262145876);

    B->b[0] = RCONST(0.740789228840);
    B->b[1] = RCONST(0.259210771159);

    B->c[0] = RCONST(0.292893218813);
    B->c[1] = RCONST(1.091883092037);
    B->c[2] = RCONST(1.292893218813);
    break;

  case(TRBDF2_3_3_2):    /* TRBDF2-ESDIRK */
    B = ARKodeButcherTable_Alloc(3, SUNTRUE);
    B->q = 2;
    B->p = 3;

    B->A[1][0] = (RCONST(2.0)-SUNRsqrt(RCONST(2.0)))/RCONST(2.0);
    B->A[1][1] = (RCONST(2.0)-SUNRsqrt(RCONST(2.0)))/RCONST(2.0);
    B->A[2][0] = SUNRsqrt(RCONST(2.0))/RCONST(4.0);
    B->A[2][1] = SUNRsqrt(RCONST(2.0))/RCONST(4.0);
    B->A[2][2] = (RCONST(2.0)-SUNRsqrt(RCONST(2.0)))/RCONST(2.0);

    B->d[0] = (RCONST(1.0)-SUNRsqrt(RCONST(2.0))/RCONST(4.0))/RCONST(3.0);
    B->d[1] = (RCONST(3.0)*SUNRsqrt(RCONST(2.0))/RCONST(4.0)+RCONST(1.0))/RCONST(3.0);
    B->d[2] = (RCONST(2.0)-SUNRsqrt(RCONST(2.0)))/RCONST(6.0);

    B->b[0] = SUNRsqrt(RCONST(2.0))/RCONST(4.0);
    B->b[1] = SUNRsqrt(RCONST(2.0))/RCONST(4.0);
    B->b[2] = (RCONST(2.0)-SUNRsqrt(RCONST(2.0)))/RCONST(2.0);

    B->c[1] = RCONST(2.0)-SUNRsqrt(RCONST(2.0));
    B->c[2] = RCONST(1.0);
    break;

  case(KVAERNO_4_2_3):    /* Kvaerno(4,2,3)-ESDIRK */
    B = ARKodeButcherTable_Alloc(4, SUNTRUE);
    B->q = 3;
    B->p = 2;
    B->A[1][0] = RCONST(0.4358665215);
    B->A[1][1] = RCONST(0.4358665215);
    B->A[2][0] = RCONST(0.490563388419108);
    B->A[2][1] = RCONST(0.073570090080892);
    B->A[2][2] = RCONST(0.4358665215);
    B->A[3][0] = RCONST(0.308809969973036);
    B->A[3][1] = RCONST(1.490563388254106);
    B->A[3][2] = RCONST(-1.235239879727145);
    B->A[3][3] = RCONST(0.4358665215);

    B->b[0] = RCONST(0.308809969973036);
    B->b[1] = RCONST(1.490563388254106);
    B->b[2] = RCONST(-1.235239879727145);
    B->b[3] = RCONST(0.4358665215);

    B->d[0] = RCONST(0.490563388419108);
    B->d[1] = RCONST(0.073570090080892);
    B->d[2] = RCONST(0.4358665215);

    B->c[1] = RCONST(0.871733043);
    B->c[2] = RCONST(1.0);
    B->c[3] = RCONST(1.0);
    break;

  case(ARK324L2SA_DIRK_4_2_3):    /* ARK3(2)4L[2]SA-ESDIRK */
    B = ARKodeButcherTable_Alloc(4, SUNTRUE);
    B->q = 3;
    B->p = 2;
    B->A[1][0] = RCONST(1767732205903.0)/RCONST(4055673282236.0);
    B->A[1][1] = RCONST(1767732205903.0)/RCONST(4055673282236.0);
    B->A[2][0] = RCONST(2746238789719.0)/RCONST(10658868560708.0);
    B->A[2][1] = RCONST(-640167445237.0)/RCONST(6845629431997.0);
    B->A[2][2] = RCONST(1767732205903.0)/RCONST(4055673282236.0);
    B->A[3][0] = RCONST(1471266399579.0)/RCONST(7840856788654.0);
    B->A[3][1] = RCONST(-4482444167858.0)/RCONST(7529755066697.0);
    B->A[3][2] = RCONST(11266239266428.0)/RCONST(11593286722821.0);
    B->A[3][3] = RCONST(1767732205903.0)/RCONST(4055673282236.0);

    B->b[0] = RCONST(1471266399579.0)/RCONST(7840856788654.0);
    B->b[1] = RCONST(-4482444167858.0)/RCONST(7529755066697.0);
    B->b[2] = RCONST(11266239266428.0)/RCONST(11593286722821.0);
    B->b[3] = RCONST(1767732205903.0)/RCONST(4055673282236.0);

    B->d[0] = RCONST(2756255671327.0)/RCONST(12835298489170.0);
    B->d[1] = RCONST(-10771552573575.0)/RCONST(22201958757719.0);
    B->d[2] = RCONST(9247589265047.0)/RCONST(10645013368117.0);
    B->d[3] = RCONST(2193209047091.0)/RCONST(5459859503100.0);

    B->c[1] = RCONST(1767732205903.0)/RCONST(2027836641118.0);
    B->c[2] = RCONST(3.0)/RCONST(5.0);
    B->c[3] = RCONST(1.0);
    break;

  case(CASH_5_2_4):    /* Cash(5,2,4)-SDIRK */
    B = ARKodeButcherTable_Alloc(5, SUNTRUE);
    B->q = 4;
    B->p = 2;
    B->A[0][0] = RCONST(0.435866521508);
    B->A[1][0] = RCONST(-1.13586652150);
    B->A[1][1] = RCONST(0.435866521508);
    B->A[2][0] = RCONST(1.08543330679);
    B->A[2][1] = RCONST(-0.721299828287);
    B->A[2][2] = RCONST(0.435866521508);
    B->A[3][0] = RCONST(0.416349501547);
    B->A[3][1] = RCONST(0.190984004184);
    B->A[3][2] = RCONST(-0.118643265417);
    B->A[3][3] = RCONST(0.435866521508);
    B->A[4][0] = RCONST(0.896869652944);
    B->A[4][1] = RCONST(0.0182725272734);
    B->A[4][2] = RCONST(-0.0845900310706);
    B->A[4][3] = RCONST(-0.266418670647);
    B->A[4][4] = RCONST(0.435866521508);

    B->b[0] = RCONST(0.896869652944);
    B->b[1] = RCONST(0.0182725272734);
    B->b[2] = RCONST(-0.0845900310706);
    B->b[3] = RCONST(-0.266418670647);
    B->b[4] = RCONST(0.435866521508);

    B->d[0] = (RCONST(-0.7)-RCONST(0.5))/(RCONST(-0.7)-RCONST(0.435866521508));
    B->d[1] = (RCONST(0.5)-RCONST(0.435866521508))/(RCONST(-0.7)-RCONST(0.435866521508));

    B->c[0] = RCONST(0.435866521508);
    B->c[1] = RCONST(-0.7);
    B->c[2] = RCONST(0.8);
    B->c[3] = RCONST(0.924556761814);
    B->c[4] = RCONST(1.0);
    break;

  case(CASH_5_3_4):    /* Cash(5,3,4)-SDIRK */
    B = ARKodeButcherTable_Alloc(5, SUNTRUE);
    B->q = 4;
    B->p = 3;
    B->A[0][0] = RCONST(0.435866521508);
    B->A[1][0] = RCONST(-1.13586652150);
    B->A[1][1] = RCONST(0.435866521508);
    B->A[2][0] = RCONST(1.08543330679);
    B->A[2][1] = RCONST(-0.721299828287);
    B->A[2][2] = RCONST(0.435866521508);
    B->A[3][0] = RCONST(0.416349501547);
    B->A[3][1] = RCONST(0.190984004184);
    B->A[3][2] = RCONST(-0.118643265417);
    B->A[3][3] = RCONST(0.435866521508);
    B->A[4][0] = RCONST(0.896869652944);
    B->A[4][1] = RCONST(0.0182725272734);
    B->A[4][2] = RCONST(-0.0845900310706);
    B->A[4][3] = RCONST(-0.266418670647);
    B->A[4][4] = RCONST(0.435866521508);

    B->b[0] = RCONST(0.896869652944);
    B->b[1] = RCONST(0.0182725272734);
    B->b[2] = RCONST(-0.0845900310706);
    B->b[3] = RCONST(-0.266418670647);
    B->b[4] = RCONST(0.435866521508);

    B->d[0] = RCONST(0.776691932910);
    B->d[1] = RCONST(0.0297472791484);
    B->d[2] = RCONST(-0.0267440239074);
    B->d[3] = RCONST(0.220304811849);

    B->c[0] = RCONST(0.435866521508);
    B->c[1] = RCONST(-0.7);
    B->c[2] = RCONST(0.8);
    B->c[3] = RCONST(0.924556761814);
    B->c[4] = RCONST(1.0);
    break;

  case(SDIRK_5_3_4):    /* SDIRK-5-4 */
    B = ARKodeButcherTable_Alloc(5, SUNTRUE);
    B->q = 4;
    B->p = 3;
    B->A[0][0] = RCONST(0.25);
    B->A[1][0] = RCONST(0.5);
    B->A[1][1] = RCONST(0.25);
    B->A[2][0] = RCONST(17.0)/RCONST(50.0);
    B->A[2][1] = RCONST(-1.0)/RCONST(25.0);
    B->A[2][2] = RCONST(0.25);
    B->A[3][0] = RCONST(371.0)/RCONST(1360.0);
    B->A[3][1] = RCONST(-137.0)/RCONST(2720.0);
    B->A[3][2] = RCONST(15.0)/RCONST(544.0);
    B->A[3][3] = RCONST(0.25);
    B->A[4][0] = RCONST(25.0)/RCONST(24.0);
    B->A[4][1] = RCONST(-49.0)/RCONST(48.0);
    B->A[4][2] = RCONST(125.0)/RCONST(16.0);
    B->A[4][3] = RCONST(-85.0)/RCONST(12.0);
    B->A[4][4] = RCONST(0.25);

    B->b[0] = RCONST(25.0)/RCONST(24.0);
    B->b[1] = RCONST(-49.0)/RCONST(48.0);
    B->b[2] = RCONST(125.0)/RCONST(16.0);
    B->b[3] = RCONST(-85.0)/RCONST(12.0);
    B->b[4] = RCONST(0.25);

    B->d[0] = RCONST(59.0)/RCONST(48.0);
    B->d[1] = RCONST(-17.0)/RCONST(96.0);
    B->d[2] = RCONST(225.0)/RCONST(32.0);
    B->d[3] = RCONST(-85.0)/RCONST(12.0);

    B->c[0] = RCONST(0.25);
    B->c[1] = RCONST(0.75);
    B->c[2] = RCONST(11.0)/RCONST(20.0);
    B->c[3] = RCONST(0.5);
    B->c[4] = RCONST(1.0);
    break;

  case(KVAERNO_5_3_4):    /* Kvaerno(5,3,4)-ESDIRK */
    B = ARKodeButcherTable_Alloc(5, SUNTRUE);
    B->q = 4;
    B->p = 3;
    B->A[1][0] = RCONST(0.4358665215);
    B->A[1][1] = RCONST(0.4358665215);
    B->A[2][0] = RCONST(0.140737774731968);
    B->A[2][1] = RCONST(-0.108365551378832);
    B->A[2][2] = RCONST(0.4358665215);
    B->A[3][0] = RCONST(0.102399400616089);
    B->A[3][1] = RCONST(-0.376878452267324);
    B->A[3][2] = RCONST(0.838612530151233);
    B->A[3][3] = RCONST(0.4358665215);
    B->A[4][0] = RCONST(0.157024897860995);
    B->A[4][1] = RCONST(0.117330441357768);
    B->A[4][2] = RCONST(0.61667803039168);
    B->A[4][3] = RCONST(-0.326899891110444);
    B->A[4][4] = RCONST(0.4358665215);

    B->b[0] = RCONST(0.157024897860995);
    B->b[1] = RCONST(0.117330441357768);
    B->b[2] = RCONST(0.61667803039168);
    B->b[3] = RCONST(-0.326899891110444);
    B->b[4] = RCONST(0.4358665215);

    B->d[0] = RCONST(0.102399400616089);
    B->d[1] = RCONST(-0.376878452267324);
    B->d[2] = RCONST(0.838612530151233);
    B->d[3] = RCONST(0.4358665215);

    B->c[1] = RCONST(0.871733043);
    B->c[2] = RCONST(0.468238744853136);
    B->c[3] = RCONST(1.0);
    B->c[4] = RCONST(1.0);
    break;

  case(ARK436L2SA_DIRK_6_3_4):    /* ARK4(3)6L[2]SA-ESDIRK */
    B = ARKodeButcherTable_Alloc(6, SUNTRUE);
    B->q = 4;
    B->p = 3;
    B->A[1][0] = RCONST(1.0)/RCONST(4.0);
    B->A[1][1] = RCONST(1.0)/RCONST(4.0);
    B->A[2][0] = RCONST(8611.0)/RCONST(62500.0);
    B->A[2][1] = RCONST(-1743.0)/RCONST(31250.0);
    B->A[2][2] = RCONST(1.0)/RCONST(4.0);
    B->A[3][0] = RCONST(5012029.0)/RCONST(34652500.0);
    B->A[3][1] = RCONST(-654441.0)/RCONST(2922500.0);
    B->A[3][2] = RCONST(174375.0)/RCONST(388108.0);
    B->A[3][3] = RCONST(1.0)/RCONST(4.0);
    B->A[4][0] = RCONST(15267082809.0)/RCONST(155376265600.0);
    B->A[4][1] = RCONST(-71443401.0)/RCONST(120774400.0);
    B->A[4][2] = RCONST(730878875.0)/RCONST(902184768.0);
    B->A[4][3] = RCONST(2285395.0)/RCONST(8070912.0);
    B->A[4][4] = RCONST(1.0)/RCONST(4.0);
    B->A[5][0] = RCONST(82889.0)/RCONST(524892.0);
    B->A[5][2] = RCONST(15625.0)/RCONST(83664.0);
    B->A[5][3] = RCONST(69875.0)/RCONST(102672.0);
    B->A[5][4] = RCONST(-2260.0)/RCONST(8211.0);
    B->A[5][5] = RCONST(1.0)/RCONST(4.0);

    B->b[0] = RCONST(82889.0)/RCONST(524892.0);
    B->b[2] = RCONST(15625.0)/RCONST(83664.0);
    B->b[3] = RCONST(69875.0)/RCONST(102672.0);
    B->b[4] = RCONST(-2260.0)/RCONST(8211.0);
    B->b[5] = RCONST(1.0)/RCONST(4.0);

    B->c[1] = RCONST(1.0)/RCONST(2.0);
    B->c[2] = RCONST(83.0)/RCONST(250.0);
    B->c[3] = RCONST(31.0)/RCONST(50.0);
    B->c[4] = RCONST(17.0)/RCONST(20.0);
    B->c[5] = RCONST(1.0);

    B->d[0] = RCONST(4586570599.0)/RCONST(29645900160.0);
    B->d[2] = RCONST(178811875.0)/RCONST(945068544.0);
    B->d[3] = RCONST(814220225.0)/RCONST(1159782912.0);
    B->d[4] = RCONST(-3700637.0)/RCONST(11593932.0);
    B->d[5] = RCONST(61727.0)/RCONST(225920.0);
    break;

  case(KVAERNO_7_4_5):    /* Kvaerno(7,4,5)-ESDIRK */
    B = ARKodeButcherTable_Alloc(7, SUNTRUE);
    B->q = 5;
    B->p = 4;
    B->A[1][0] = RCONST(0.26);
    B->A[1][1] = RCONST(0.26);
    B->A[2][0] = RCONST(0.13);
    B->A[2][1] = RCONST(0.84033320996790809);
    B->A[2][2] = RCONST(0.26);
    B->A[3][0] = RCONST(0.22371961478320505);
    B->A[3][1] = RCONST(0.47675532319799699);
    B->A[3][2] = RCONST(-0.06470895363112615);
    B->A[3][3] = RCONST(0.26);
    B->A[4][0] = RCONST(0.16648564323248321);
    B->A[4][1] = RCONST(0.10450018841591720);
    B->A[4][2] = RCONST(0.03631482272098715);
    B->A[4][3] = RCONST(-0.13090704451073998);
    B->A[4][4] = RCONST(0.26);
    B->A[5][0] = RCONST(0.13855640231268224);
    B->A[5][2] = RCONST(-0.04245337201752043);
    B->A[5][3] = RCONST(0.02446657898003141);
    B->A[5][4] = RCONST(0.61943039072480676);
    B->A[5][5] = RCONST(0.26);
    B->A[6][0] = RCONST(0.13659751177640291);
    B->A[6][2] = RCONST(-0.05496908796538376);
    B->A[6][3] = RCONST(-0.04118626728321046);
    B->A[6][4] = RCONST(0.62993304899016403);
    B->A[6][5] = RCONST(0.06962479448202728);
    B->A[6][6] = RCONST(0.26);

    B->b[0] = RCONST(0.13659751177640291);
    B->b[2] = RCONST(-0.05496908796538376);
    B->b[3] = RCONST(-0.04118626728321046);
    B->b[4] = RCONST(0.62993304899016403);
    B->b[5] = RCONST(0.06962479448202728);
    B->b[6] = RCONST(0.26);

    B->d[0] = RCONST(0.13855640231268224);
    B->d[2] = RCONST(-0.04245337201752043);
    B->d[3] = RCONST(0.02446657898003141);
    B->d[4] = RCONST(0.61943039072480676);
    B->d[5] = RCONST(0.26);

    B->c[1] = RCONST(0.52);
    B->c[2] = RCONST(1.230333209967908);
    B->c[3] = RCONST(0.895765984350076);
    B->c[4] = RCONST(0.436393609858648);
    B->c[5] = RCONST(1.0);
    B->c[6] = RCONST(1.0);
    break;

  case(ARK548L2SA_DIRK_8_4_5):    /* ARK5(4)8L[2]SA-ESDIRK */
    B = ARKodeButcherTable_Alloc(8, SUNTRUE);
    B->q = 5;
    B->p = 4;
    B->A[1][0] = RCONST(41.0)/RCONST(200.0);
    B->A[1][1] = RCONST(41.0)/RCONST(200.0);
    B->A[2][0] = RCONST(41.0)/RCONST(400.0);
    B->A[2][1] = RCONST(-567603406766.0)/RCONST(11931857230679.0);
    B->A[2][2] = RCONST(41.0)/RCONST(200.0);
    B->A[3][0] = RCONST(683785636431.0)/RCONST(9252920307686.0);
    B->A[3][2] = RCONST(-110385047103.0)/RCONST(1367015193373.0);
    B->A[3][3] = RCONST(41.0)/RCONST(200.0);
    B->A[4][0] = RCONST(3016520224154.0)/RCONST(10081342136671.0);
    B->A[4][2] = RCONST(30586259806659.0)/RCONST(12414158314087.0);
    B->A[4][3] = RCONST(-22760509404356.0)/RCONST(11113319521817.0);
    B->A[4][4] = RCONST(41.0)/RCONST(200.0);
    B->A[5][0] = RCONST(218866479029.0)/RCONST(1489978393911.0);
    B->A[5][2] = RCONST(638256894668.0)/RCONST(5436446318841.0);
    B->A[5][3] = RCONST(-1179710474555.0)/RCONST(5321154724896.0);
    B->A[5][4] = RCONST(-60928119172.0)/RCONST(8023461067671.0);
    B->A[5][5] = RCONST(41.0)/RCONST(200.0);
    B->A[6][0] = RCONST(1020004230633.0)/RCONST(5715676835656.0);
    B->A[6][2] = RCONST(25762820946817.0)/RCONST(25263940353407.0);
    B->A[6][3] = RCONST(-2161375909145.0)/RCONST(9755907335909.0);
    B->A[6][4] = RCONST(-211217309593.0)/RCONST(5846859502534.0);
    B->A[6][5] = RCONST(-4269925059573.0)/RCONST(7827059040749.0);
    B->A[6][6] = RCONST(41.0)/RCONST(200.0);
    B->A[7][0] = RCONST(-872700587467.0)/RCONST(9133579230613.0);
    B->A[7][3] = RCONST(22348218063261.0)/RCONST(9555858737531.0);
    B->A[7][4] = RCONST(-1143369518992.0)/RCONST(8141816002931.0);
    B->A[7][5] = RCONST(-39379526789629.0)/RCONST(19018526304540.0);
    B->A[7][6] = RCONST(32727382324388.0)/RCONST(42900044865799.0);
    B->A[7][7] = RCONST(41.0)/RCONST(200.0);

    B->b[0] = RCONST(-872700587467.0)/RCONST(9133579230613.0);
    B->b[3] = RCONST(22348218063261.0)/RCONST(9555858737531.0);
    B->b[4] = RCONST(-1143369518992.0)/RCONST(8141816002931.0);
    B->b[5] = RCONST(-39379526789629.0)/RCONST(19018526304540.0);
    B->b[6] = RCONST(32727382324388.0)/RCONST(42900044865799.0);
    B->b[7] = RCONST(41.0)/RCONST(200.0);

    B->d[0] = RCONST(-975461918565.0)/RCONST(9796059967033.0);
    B->d[3] = RCONST(78070527104295.0)/RCONST(32432590147079.0);
    B->d[4] = RCONST(-548382580838.0)/RCONST(3424219808633.0);
    B->d[5] = RCONST(-33438840321285.0)/RCONST(15594753105479.0);
    B->d[6] = RCONST(3629800801594.0)/RCONST(4656183773603.0);
    B->d[7] = RCONST(4035322873751.0)/RCONST(18575991585200.0);

    B->c[1] = RCONST(41.0)/RCONST(100.0);
    B->c[2] = RCONST(2935347310677.0)/RCONST(11292855782101.0);
    B->c[3] = RCONST(1426016391358.0)/RCONST(7196633302097.0);
    B->c[4] = RCONST(92.0)/RCONST(100.0);
    B->c[5] = RCONST(24.0)/RCONST(100.0);
    B->c[6] = RCONST(3.0)/RCONST(5.0);
    B->c[7] = RCONST(1.0);
    break;

  default:

    arkProcessError(NULL, ARK_ILL_INPUT, "ARKode",
                    "ARKodeButcherTable_LoadDIRK",
                    "Unknown Butcher table");
    return(NULL);

  }

  return(B);
}


/*---------------------------------------------------------------
  EOF
  ---------------------------------------------------------------*/
