/* ---------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 *                Daniel R. Reynolds @ SMU
 * ---------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * ---------------------------------------------------------------------------
 * This is the implementation file for ARKODE's MRIStepCoupling tables.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "arkode_impl.h"
#include "arkode_mristep_impl.h"
#include <sundials/sundials_math.h>


/* ===========================================================================
 * Exported Functions
 * ===========================================================================*/


/*---------------------------------------------------------------
  Returns MRIStepCoupling table structure for pre-set MRI methods.

  Input:  imeth -- integer key for the desired method (see below)

  Allowed 'method' names and properties are listed in the table
  below.

  The 'type' column denotes whether the method is explicit (E),
  or solve-decoupled implicit (ID).

  The 'QP' column denotes whether the coefficients of the method
  are known precisely enough for use in quad precision (128-bit)
  calculations.

     imeth                       order   type    QP
    ------------------------------------------------
     ARKODE_MIS_KW3                     3       E       Y
     ARKODE_MRI_GARK_ERK33a             3       E       Y
     ARKODE_MRI_GARK_ERK45a             4       E       Y
     ARKODE_MRI_GARK_IRK21a             2       ID      Y
     ARKODE_MRI_GARK_ESDIRK34a          3       ID      Y
     ARKODE_MRI_GARK_ESDIRK46a          4       ID      Y
     ARKODE_IMEX_MRI_GARK3a             3       ID      Y
     ARKODE_IMEX_MRI_GARK3b             3       ID      Y
     ARKODE_IMEX_MRI_GARK4              4       ID      Y
    ------------------------------------------------

  ---------------------------------------------------------------*/
MRIStepCoupling MRIStepCoupling_LoadTable(ARKODE_MRITableID imethod)
{

  MRIStepCoupling    C = NULL;
  ARKodeButcherTable B = NULL;
  realtype beta;

  /* fill in coefficients based on method name */
  switch(imethod) {

  case(ARKODE_MIS_KW3):
    /* Schlegel et al., JCAM 226:345-357, 2009 */
    B = ARKodeButcherTable_LoadERK(ARKODE_KNOTH_WOLKE_3_3);
    C = MRIStepCoupling_MIStoMRI(B, 3, 0);
    ARKodeButcherTable_Free(B);
    break;

  case(ARKODE_MRI_GARK_ERK33a):
    /* A. Sandu, SINUM 57:2300-2327, 2019 */
    C = MRIStepCoupling_Alloc(2, 4, MRISTEP_EXPLICIT);

    C->q = 3;
    C->p = 0;

    C->c[1] = ONE/RCONST(3.0);
    C->c[2] = TWO/RCONST(3.0);
    C->c[3] = ONE;

    C->W[0][1][0] =  ONE/RCONST(3.0);
    C->W[0][2][0] = -ONE/RCONST(3.0);
    C->W[0][2][1] =  TWO/RCONST(3.0);
    C->W[0][3][1] = -TWO/RCONST(3.0);
    C->W[0][3][2] =  ONE;

    C->W[1][3][0] =  ONE/TWO;
    C->W[1][3][2] = -ONE/TWO;

    break;

  case(ARKODE_MRI_GARK_ERK45a):
    /* A. Sandu, SINUM 57:2300-2327, 2019 */
    C = MRIStepCoupling_Alloc(2, 6, MRISTEP_EXPLICIT);

    C->q = 4;
    C->p = 0;

    C->c[1] = RCONST(0.2);
    C->c[2] = RCONST(0.4);
    C->c[3] = RCONST(0.6);
    C->c[4] = RCONST(0.8);
    C->c[5] = ONE;

    C->W[0][1][0] =  RCONST(0.2);
    C->W[0][2][0] = -RCONST(53.0)/RCONST(16.0);
    C->W[0][2][1] =  RCONST(281.0)/RCONST(80.0);
    C->W[0][3][0] = -RCONST(36562993.0)/RCONST(71394880.0);
    C->W[0][3][1] =  RCONST(34903117.0)/RCONST(17848720.0);
    C->W[0][3][2] = -RCONST(88770499.0)/RCONST(71394880.0);
    C->W[0][4][0] = -RCONST(7631593.0)/RCONST(71394880.0);
    C->W[0][4][1] = -RCONST(166232021.0)/RCONST(35697440.0);
    C->W[0][4][2] =  RCONST(6068517.0)/RCONST(1519040.0);
    C->W[0][4][3] =  RCONST(8644289.0)/RCONST(8924360.0);
    C->W[0][5][0] =  RCONST(277061.0)/RCONST(303808.0);
    C->W[0][5][1] = -RCONST(209323.0)/RCONST(1139280.0);
    C->W[0][5][2] = -RCONST(1360217.0)/RCONST(1139280.0);
    C->W[0][5][3] = -RCONST(148789.0)/RCONST(56964.0);
    C->W[0][5][4] =  RCONST(147889.0)/RCONST(45120.0);

    C->W[1][2][0] =  RCONST(503.0)/RCONST(80.0);
    C->W[1][2][1] = -RCONST(503.0)/RCONST(80.0);
    C->W[1][3][0] = -RCONST(1365537.0)/RCONST(35697440.0);
    C->W[1][3][1] =  RCONST(4963773.0)/RCONST(7139488.0);
    C->W[1][3][2] = -RCONST(1465833.0)/RCONST(2231090.0);
    C->W[1][4][0] =  RCONST(66974357.0)/RCONST(35697440.0);
    C->W[1][4][1] =  RCONST(21445367.0)/RCONST(7139488.0);
    C->W[1][4][2] = -RCONST(3.0);
    C->W[1][4][3] = -RCONST(8388609.0)/RCONST(4462180.0);
    C->W[1][5][0] = -RCONST(18227.0)/RCONST(7520.0);
    C->W[1][5][1] =  TWO;
    C->W[1][5][2] =  ONE;
    C->W[1][5][3] =  RCONST(5.0);
    C->W[1][5][4] = -RCONST(41933.0)/RCONST(7520.0);

    break;

  case(ARKODE_MRI_GARK_IRK21a):
    /* A. Sandu, SINUM 57:2300-2327, 2019 */
    B = ARKodeButcherTable_Alloc(3, SUNFALSE);

    B->q=2;

    B->c[1] = ONE;
    B->c[2] = ONE;

    B->A[1][0] = ONE;
    B->A[2][0] = RCONST(0.5);
    B->A[2][2] = RCONST(0.5);

    B->b[0] = RCONST(0.5);
    B->b[2] = RCONST(0.5);

    C = MRIStepCoupling_MIStoMRI(B, 2, 0);
    ARKodeButcherTable_Free(B);

    break;

  case(ARKODE_MRI_GARK_ESDIRK34a):
    /* A. Sandu, SINUM 57:2300-2327, 2019 */
    C = MRIStepCoupling_Alloc(1, 7, MRISTEP_IMPLICIT);

    beta = RCONST(0.4358665215084589994160194511935568425);

    C->q = 3;
    C->p = 0;

    C->c[1] = ONE/RCONST(3.0);
    C->c[2] = ONE/RCONST(3.0);
    C->c[3] = TWO/RCONST(3.0);
    C->c[4] = TWO/RCONST(3.0);
    C->c[5] = ONE;
    C->c[6] = ONE;

    C->G[0][1][0] =  ONE/RCONST(3.0);
    C->G[0][2][0] = -beta;
    C->G[0][2][2] =  beta;
    C->G[0][3][0] =  RCONST(-0.3045790611944504970424837655380884888);
    C->G[0][3][2] =  RCONST(0.6379123945277838303758170988714218222);
    C->G[0][4][0] =  RCONST(0.2116913105640266601676536489364004869);
    C->G[0][4][2] =  RCONST(-0.6475578320724856595836731001299573294);
    C->G[0][4][4] =  beta;
    C->G[0][5][0] =  RCONST(0.4454209388055495029575162344619115112);
    C->G[0][5][2] =  RCONST(0.8813784805616198280398949036456491923);
    C->G[0][5][4] =  RCONST(-0.9934660860338359976640778047742273701);
    C->G[0][6][0] = -beta;
    C->G[0][6][6] =  beta;

    break;

  case(ARKODE_MRI_GARK_ESDIRK46a):
    /* A. Sandu, SINUM 57:2300-2327, 2019 */
    C = MRIStepCoupling_Alloc(2, 11, MRISTEP_IMPLICIT);

    C->q = 4;
    C->p = 0;

    C->c[1]  = ONE/RCONST(5.0);
    C->c[2]  = ONE/RCONST(5.0);
    C->c[3]  = TWO/RCONST(5.0);
    C->c[4]  = TWO/RCONST(5.0);
    C->c[5]  = RCONST(3.0)/RCONST(5.0);
    C->c[6]  = RCONST(3.0)/RCONST(5.0);
    C->c[7]  = RCONST(4.0)/RCONST(5.0);
    C->c[8]  = RCONST(4.0)/RCONST(5.0);
    C->c[9]  = ONE;
    C->c[10] = ONE;

    C->G[0][1][0]   =  ONE/RCONST(5.0);
    C->G[0][2][0]   = -ONE/RCONST(4.0);
    C->G[0][2][2]   =  ONE/RCONST(4.0);
    C->G[0][3][0]   =  RCONST(1771023115159.0)/RCONST(1929363690800.0);
    C->G[0][3][2]   = -RCONST(1385150376999.0)/RCONST(1929363690800.0);
    C->G[0][4][0]   =  RCONST(914009.0)/RCONST(345800.0);
    C->G[0][4][2]   = -RCONST(1000459.0)/RCONST(345800.0);
    C->G[0][4][4]   =  ONE/RCONST(4.0);
    C->G[0][5][0]   =  RCONST(18386293581909.0)/RCONST(36657910125200.0);
    C->G[0][5][2]   =  RCONST(5506531089.0)/RCONST(80566835440.0);
    C->G[0][5][4]   = -RCONST(178423463189.0)/RCONST(482340922700.0);
    C->G[0][6][0]   =  RCONST(36036097.0)/RCONST(8299200.0);
    C->G[0][6][2]   =  RCONST(4621.0)/RCONST(118560.0);
    C->G[0][6][4]   = -RCONST(38434367.0)/RCONST(8299200.0);
    C->G[0][6][6]   =  ONE/RCONST(4.0);
    C->G[0][7][0]   = -RCONST(247809665162987.0)/RCONST(146631640500800.0);
    C->G[0][7][2]   =  RCONST(10604946373579.0)/RCONST(14663164050080.0);
    C->G[0][7][4]   =  RCONST(10838126175385.0)/RCONST(5865265620032.0);
    C->G[0][7][6]   = -RCONST(24966656214317.0)/RCONST(36657910125200.0);
    C->G[0][8][0]   =  RCONST(38519701.0)/RCONST(11618880.0);
    C->G[0][8][2]   =  RCONST(10517363.0)/RCONST(9682400.0);
    C->G[0][8][4]   = -RCONST(23284701.0)/RCONST(19364800.0);
    C->G[0][8][6]   = -RCONST(10018609.0)/RCONST(2904720.0);
    C->G[0][8][8]   =  ONE/RCONST(4.0);
    C->G[0][9][0]   = -RCONST(52907807977903.0)/RCONST(33838070884800.0);
    C->G[0][9][2]   =  RCONST(74846944529257.0)/RCONST(73315820250400.0);
    C->G[0][9][4]   =  RCONST(365022522318171.0)/RCONST(146631640500800.0);
    C->G[0][9][6]   = -RCONST(20513210406809.0)/RCONST(109973730375600.0);
    C->G[0][9][8]   = -RCONST(2918009798.0)/RCONST(1870301537.0);
    C->G[0][10][0]  =  RCONST(19.0)/RCONST(100.0);
    C->G[0][10][2]  = -RCONST(73.0)/RCONST(300.0);
    C->G[0][10][4]  =  RCONST(127.0)/RCONST(300.0);
    C->G[0][10][6]  =  RCONST(127.0)/RCONST(300.0);
    C->G[0][10][8]  = -RCONST(313.0)/RCONST(300.0);
    C->G[0][10][10] =  ONE/RCONST(4.0);

    C->G[1][3][0]  = -RCONST(1674554930619.0)/RCONST(964681845400.0);
    C->G[1][3][2]  =  RCONST(1674554930619.0)/RCONST(964681845400.0);
    C->G[1][4][0]  = -RCONST(1007739.0)/RCONST(172900.0);
    C->G[1][4][2]  =  RCONST(1007739.0)/RCONST(172900.0);
    C->G[1][5][0]  = -RCONST(8450070574289.0)/RCONST(18328955062600.0);
    C->G[1][5][2]  = -RCONST(39429409169.0)/RCONST(40283417720.0);
    C->G[1][5][4]  =  RCONST(173621393067.0)/RCONST(120585230675.0);
    C->G[1][6][0]  = -RCONST(122894383.0)/RCONST(16598400.0);
    C->G[1][6][2]  =  RCONST(14501.0)/RCONST(237120.0);
    C->G[1][6][4]  =  RCONST(121879313.0)/RCONST(16598400.0);
    C->G[1][7][0]  =  RCONST(32410002731287.0)/RCONST(15434909526400.0);
    C->G[1][7][2]  = -RCONST(46499276605921.0)/RCONST(29326328100160.0);
    C->G[1][7][4]  = -RCONST(34914135774643.0)/RCONST(11730531240064.0);
    C->G[1][7][6]  =  RCONST(45128506783177.0)/RCONST(18328955062600.0);
    C->G[1][8][0]  = -RCONST(128357303.0)/RCONST(23237760.0);
    C->G[1][8][2]  = -RCONST(35433927.0)/RCONST(19364800.0);
    C->G[1][8][4]  =  RCONST(71038479.0)/RCONST(38729600.0);
    C->G[1][8][6]  =  RCONST(8015933.0)/RCONST(1452360.0);
    C->G[1][9][0]  =  RCONST(136721604296777.0)/RCONST(67676141769600.0);
    C->G[1][9][2]  = -RCONST(349632444539303.0)/RCONST(146631640500800.0);
    C->G[1][9][4]  = -RCONST(1292744859249609.0)/RCONST(293263281001600.0);
    C->G[1][9][6]  =  RCONST(8356250416309.0)/RCONST(54986865187800.0);
    C->G[1][9][8]  =  RCONST(17282943803.0)/RCONST(3740603074.0);
    C->G[1][10][0] =  RCONST(3.0)/RCONST(25.0);
    C->G[1][10][2] = -RCONST(29.0)/RCONST(300.0);
    C->G[1][10][4] =  RCONST(71.0)/RCONST(300.0);
    C->G[1][10][6] =  RCONST(71.0)/RCONST(300.0);
    C->G[1][10][8] = -RCONST(149.0)/RCONST(300.0);

    break;

  case(ARKODE_IMEX_MRI_GARK3a):
    /* R. Chinomona & D. Reynolds SINUM 43(5):A3082-A3113, 2021 */
    C = MRIStepCoupling_Alloc(1, 8, MRISTEP_IMEX);

    beta = RCONST(0.4358665215084589994160194511935568425);

    C->q = 3;
    C->p = 0;

    C->c[1] = beta;
    C->c[2] = beta;
    C->c[3] = RCONST(0.7179332607542294997080097255967784213);
    C->c[4] = RCONST(0.7179332607542294997080097255967784213);
    C->c[5] = ONE;
    C->c[6] = ONE;
    C->c[7] = ONE;

    C->W[0][1][0] =  beta;
    C->W[0][3][0] = -RCONST(0.5688715801234400928465032925317932021);
    C->W[0][3][2] =  RCONST(0.8509383193692105931384935669350147809);
    C->W[0][4][0] =  RCONST(0.454283944643608855878770886900124654);
    C->W[0][4][2] = -RCONST(0.454283944643608855878770886900124654);
    C->W[0][5][0] = -RCONST(0.4271371821005074011706645050390732474);
    C->W[0][5][2] =  RCONST(0.1562747733103380821014660497037023496);
    C->W[0][5][4] =  RCONST(0.5529291480359398193611887297385924765);
    C->W[0][7][0] =  RCONST(0.105858296071879638722377459477184953);
    C->W[0][7][2] =  RCONST(0.655567501140070250975288954324730635);
    C->W[0][7][4] = -RCONST(1.197292318720408889113685864995472431);
    C->W[0][7][6] =  beta;

    C->G[0][1][0] =  beta;
    C->G[0][2][0] = -beta;
    C->G[0][2][2] =  beta;
    C->G[0][3][0] = -RCONST(0.4103336962288525014599513720161078937);
    C->G[0][3][2] =  RCONST(0.6924004354746230017519416464193294724);
    C->G[0][4][0] =  RCONST(0.4103336962288525014599513720161078937);
    C->G[0][4][2] = -RCONST(0.8462002177373115008759708232096647362);
    C->G[0][4][4] =  beta;
    C->G[0][5][0] =  beta;
    C->G[0][5][2] =  RCONST(0.9264299099302395700444874096601015328);
    C->G[0][5][4] = -RCONST(1.080229692192928069168516586450436797);
    C->G[0][6][0] = -beta;
    C->G[0][6][6] =  beta;

    break;

  case(ARKODE_IMEX_MRI_GARK3b):
    /* R. Chinomona & D. Reynolds SINUM 43(5):A3082-A3113, 2021 */
    C = MRIStepCoupling_Alloc(1, 8, MRISTEP_IMEX);

    beta = RCONST(0.4358665215084589994160194511935568425);

    C->q = 3;
    C->p = 0;

    C->c[1] = beta;
    C->c[2] = beta;
    C->c[3] = RCONST(0.7179332607542294997080097255967784213);
    C->c[4] = RCONST(0.7179332607542294997080097255967784213);
    C->c[5] = ONE;
    C->c[6] = ONE;
    C->c[7] = ONE;

    C->W[0][1][0] =  beta;
    C->W[0][3][0] = -RCONST(0.1750145285570467590610670000018749059);
    C->W[0][3][2] =  RCONST(0.4570812678028172593530572744050964846);
    C->W[0][4][0] =  RCONST(0.06042689307721552209333459437020635774);
    C->W[0][4][2] = -RCONST(0.06042689307721552209333459437020635774);
    C->W[0][5][0] =  RCONST(0.1195213959425454440038786034027936869);
    C->W[0][5][2] = -RCONST(1.84372522668966191789853395029629765);
    C->W[0][5][4] =  RCONST(2.006270569992886974186645621296725542);
    C->W[0][6][0] = -RCONST(0.5466585780430528451745431084418669343);
    C->W[0][6][2] =  RCONST(2.0);
    C->W[0][6][4] = -RCONST(1.453341421956947154825456891558133066);
    C->W[0][7][0] =  RCONST(0.105858296071879638722377459477184953);
    C->W[0][7][2] =  RCONST(0.655567501140070250975288954324730635);
    C->W[0][7][4] = -RCONST(1.197292318720408889113685864995472431);
    C->W[0][7][6] =  beta;

    C->G[0][1][0] =  beta;
    C->G[0][2][0] = -beta;
    C->G[0][2][2] =  beta;
    C->G[0][3][0] =  RCONST(0.0414273753564414837153799230278275639);
    C->G[0][3][2] =  RCONST(0.2406393638893290165766103513753940148);
    C->G[0][4][0] = -RCONST(0.0414273753564414837153799230278275639);
    C->G[0][4][2] = -RCONST(0.3944391461520175157006395281657292786);
    C->G[0][4][4] =  beta;
    C->G[0][5][0] =  RCONST(0.1123373143006047802633543416889605123);
    C->G[0][5][2] =  RCONST(1.051807513648115027700693049638099167);
    C->G[0][5][4] = -RCONST(0.8820780887029493076720571169238381009);
    C->G[0][6][0] = -RCONST(0.1123373143006047802633543416889605123);
    C->G[0][6][2] = -RCONST(0.1253776037178754576562056399779976346);
    C->G[0][6][4] = -RCONST(0.1981516034899787614964594695265986957);
    C->G[0][6][6] =  beta;

    break;

  case(ARKODE_IMEX_MRI_GARK4):
    /* R. Chinomona & D. Reynolds SINUM 43(5):A3082-A3113, 2021 */
    C = MRIStepCoupling_Alloc(2, 12, MRISTEP_IMEX);

    C->q = 4;
    C->p = 0;

    C->c[1]  = RCONST(0.5);
    C->c[2]  = RCONST(0.5);
    C->c[3]  = RCONST(0.625);
    C->c[4]  = RCONST(0.625);
    C->c[5]  = RCONST(0.75);
    C->c[6]  = RCONST(0.75);
    C->c[7]  = RCONST(0.875);
    C->c[8]  = RCONST(0.875);
    C->c[9]  = ONE;
    C->c[10] = ONE;
    C->c[11] = ONE;

    C->W[0][1][0]   =  RCONST(0.5);
    C->W[0][3][0]   = -RCONST(1.91716534363662868878172216064946905);
    C->W[0][3][2]   =  RCONST(2.04216534363662868878172216064946905);
    C->W[0][4][0]   = -RCONST(0.4047510318011059426979159070469904691);
    C->W[0][4][2]   =  RCONST(0.4047510318011059426979159070469904691);
    C->W[0][5][0]   =  RCONST(11.45146602249221636665698028602631728);
    C->W[0][5][2]   = -RCONST(30.21075747526504271440647815573950607);
    C->W[0][5][4]   =  RCONST(18.88429145277282634774949786971318879);
    C->W[0][6][0]   = -RCONST(0.7090335647602614506847116729463301439);
    C->W[0][6][2]   =  RCONST(1.03030720858751876652616190884004718);
    C->W[0][6][4]   = -RCONST(0.3212736438272573158414502358937170357);
    C->W[0][7][0]   = -RCONST(29.99548716455828439840910684944199275);
    C->W[0][7][2]   =  RCONST(37.60598277499180180536489685624385701);
    C->W[0][7][4]   =  RCONST(0.3212736438272573158414502358937170357);
    C->W[0][7][6]   = -RCONST(7.806769254260774722797240242695581295);
    C->W[0][8][0]   =  RCONST(3.104665054272962116338769391849124223);
    C->W[0][8][2]   = -RCONST(2.430325019757162297132065927415566359);
    C->W[0][8][4]   = -RCONST(1.905479301151524635219201659483842131);
    C->W[0][8][6]   =  RCONST(1.231139266635724816012498195050284266);
    C->W[0][9][0]   = -RCONST(2.424429547752047869875875914355514008);
    C->W[0][9][2]   =  RCONST(2.430325019757162297132065927415566359);
    C->W[0][9][4]   =  RCONST(1.905479301151524635219201659483842131);
    C->W[0][9][6]   = -RCONST(1.231139266635724816012498195050284266);
    C->W[0][9][8]   = -RCONST(0.555235506520914246462893477493610215);
    C->W[0][10][0]  = -RCONST(0.01044135044479748590294518945165354204);
    C->W[0][10][2]  =  RCONST(0.07260303614655074505152104505488141613);
    C->W[0][10][4]  = -RCONST(0.1288275951677260952239454098576424313);
    C->W[0][10][6]  =  RCONST(0.1129355350093823566139440107122154084);
    C->W[0][10][8]  = -RCONST(0.04626962554340952053857445645780085125);
    C->W[0][11][0]  = -RCONST(0.8108522787762101328175789228607932098);
    C->W[0][11][2]  =  RCONST(0.2560073199220492435001562192140882299);
    C->W[0][11][4]  =  RCONST(0.8068294072697527893665866422787819475);
    C->W[0][11][6]  = -RCONST(0.4557148228721823795105894821742761164);
    C->W[0][11][8]  = -RCONST(0.04626962554340952053857445645780085125);
    C->W[0][11][10] =  RCONST(0.25);

    C->W[1][3][0]   =  RCONST(4.084330687273257377563444321298938099);
    C->W[1][3][2]   = -RCONST(4.084330687273257377563444321298938099);
    C->W[1][5][0]   = -RCONST(21.84342998138222084791812875795865363);
    C->W[1][5][2]   =  RCONST(59.61201288692787354341712449738503121);
    C->W[1][5][4]   = -RCONST(37.76858290554565269549899573942637758);
    C->W[1][7][0]   =  RCONST(61.65904145863709169818763704477664579);
    C->W[1][7][2]   = -RCONST(77.27257996715864114378211753016780838);
    C->W[1][7][6]   =  RCONST(15.61353850852154944559448048539116259);
    C->W[1][9][0]   = -RCONST(1.11047101304182849292578695498722043);
    C->W[1][9][8]   =  RCONST(1.11047101304182849292578695498722043);

    C->G[0][1][0]   =  RCONST(0.5);
    C->G[0][2][0]   = -RCONST(0.25);
    C->G[0][2][2]   =  RCONST(0.25);
    C->G[0][3][0]   = -RCONST(3.977281248108488183067033851462278892);
    C->G[0][3][2]   =  RCONST(4.102281248108488183067033851462278892);
    C->G[0][4][0]   = -RCONST(0.06905388741401691232724147084809374064);
    C->G[0][4][2]   = -RCONST(0.1809461125859830876727585291519062594);
    C->G[0][4][4]   =  RCONST(0.25);
    C->G[0][5][0]   = -RCONST(1.761767663757920528863378964822412405);
    C->G[0][5][2]   =  RCONST(2.694524698377298610155338150791461384);
    C->G[0][5][4]   = -RCONST(0.8077570346193780812919591859690489783);
    C->G[0][6][0]   =  RCONST(0.5558721791553969487305081009588084962);
    C->G[0][6][2]   = -RCONST(0.6799140501579995013958501527883486949);
    C->G[0][6][4]   = -RCONST(0.1259581289973974473346579481704598013);
    C->G[0][6][6]   =  RCONST(0.25);
    C->G[0][7][0]   = -RCONST(5.840176028724955954446426657541065113);
    C->G[0][7][2]   =  RCONST(8.174456684291915089191270805710716374);
    C->G[0][7][4]   =  RCONST(0.1259581289973974473346579481704598013);
    C->G[0][7][6]   = -RCONST(2.335238784564356582079502096340111063);
    C->G[0][8][0]   = -RCONST(1.906792645167811808094759305036052304);
    C->G[0][8][2]   = -RCONST(1.547057811385123933632984579249388443);
    C->G[0][8][4]   =  RCONST(4.129888013149350305954491738020313225);
    C->G[0][8][6]   = -RCONST(0.9260375565964145642267478537348724775);
    C->G[0][8][8]   =  RCONST(0.25);
    C->G[0][9][0]   =  RCONST(3.337028151688726054557652782529662519);
    C->G[0][9][2]   =  RCONST(1.547057811385123933632984579249388443);
    C->G[0][9][4]   = -RCONST(4.129888013149350305954491738020313225);
    C->G[0][9][6]   =  RCONST(0.9260375565964145642267478537348724775);
    C->G[0][9][8]   = -RCONST(1.555235506520914246462893477493610215);
    C->G[0][10][0]  = -RCONST(0.8212936292210076187205241123124467518);
    C->G[0][10][2]  =  RCONST(0.328610356068599988551677264268969646);
    C->G[0][10][4]  =  RCONST(0.6780018121020266941426412324211395162);
    C->G[0][10][6]  = -RCONST(0.3427792878628000228966454714620607079);
    C->G[0][10][8]  = -RCONST(0.0925392510868190410771489129156017025);
    C->G[0][10][10] =  RCONST(0.25);

    C->G[1][3][0]   =  RCONST(8.704562496216976366134067702924557783);
    C->G[1][3][2]   = -RCONST(8.704562496216976366134067702924557783);
    C->G[1][5][0]   =  RCONST(3.911643102343874882381240871341012292);
    C->G[1][5][2]   = -RCONST(5.027157171582631044965159243279110249);
    C->G[1][5][4]   =  RCONST(1.115514069238756162583918371938097957);
    C->G[1][7][0]   =  RCONST(10.81860769913911801143183711316451323);
    C->G[1][7][2]   = -RCONST(14.98908526826783117559084130584473536);
    C->G[1][7][6]   =  RCONST(4.170477569128713164159004192680222125);
    C->G[1][9][0]   = -RCONST(2.61047101304182849292578695498722043);
    C->G[1][9][8]   =  RCONST(2.61047101304182849292578695498722043);

    break;

  default:
    arkProcessError(NULL, ARK_ILL_INPUT, "ARKode",
                    "MRIStepCoupling_LoadTable",
                    "Unknown coupling table");
    return(NULL);
  }

  return(C);
}


/*---------------------------------------------------------------
  Routine to allocate an empty MRIStepCoupling structure
  ---------------------------------------------------------------*/
MRIStepCoupling MRIStepCoupling_Alloc(int nmat, int stages,
                                      MRISTEP_METHOD_TYPE type)
{
  int i, j;
  MRIStepCoupling MRIC = NULL;

  /* Check for legal input values */
  if (nmat < 1 || stages < 1) return(NULL);

  /* ------------------------------------------
   * Allocate and initialize coupling structure
   * ------------------------------------------ */

  MRIC = (MRIStepCoupling) malloc(sizeof(struct MRIStepCouplingMem));
  if (!MRIC) return(NULL);

  MRIC->nmat   = nmat;
  MRIC->stages = stages;
  MRIC->q      = 0;
  MRIC->p      = 0;
  MRIC->c      = NULL;
  MRIC->W      = NULL;
  MRIC->G      = NULL;

  /* --------------------------------------------
   * Allocate abscissae and coupling coefficients
   * -------------------------------------------- */

  MRIC->c = (realtype *) calloc( stages, sizeof(realtype) );
  if (!(MRIC->c)) { MRIStepCoupling_Free(MRIC); return(NULL); }

  if (type == MRISTEP_EXPLICIT || type == MRISTEP_IMEX) {

    /* allocate W matrices */
    MRIC->W = (realtype ***) calloc( nmat, sizeof(realtype**) );
    if (!(MRIC->W)) { MRIStepCoupling_Free(MRIC); return(NULL); }

    /* allocate rows of each matrix in W */
    for (i=0; i<nmat; i++) {
      MRIC->W[i] = NULL;
      MRIC->W[i] = (realtype **) calloc( stages, sizeof(realtype*) );
      if (!(MRIC->W[i])) { MRIStepCoupling_Free(MRIC); return(NULL); }
    }

    /* allocate columns of each matrix in W */
    for (i=0; i<nmat; i++)
      for (j=0; j<stages; j++) {
        MRIC->W[i][j] = NULL;
        MRIC->W[i][j] = (realtype *) calloc( stages, sizeof(realtype) );
        if (!(MRIC->W[i][j])) { MRIStepCoupling_Free(MRIC); return(NULL); }
      }
  }

  if (type == MRISTEP_IMPLICIT || type == MRISTEP_IMEX) {

    /* allocate G matrices */
    MRIC->G = (realtype ***) calloc( nmat, sizeof(realtype**) );
    if (!(MRIC->G)) { MRIStepCoupling_Free(MRIC); return(NULL); }

    /* allocate rows of each matrix in G */
    for (i=0; i<nmat; i++) {
      MRIC->G[i] = NULL;
      MRIC->G[i] = (realtype **) calloc( stages, sizeof(realtype*) );
      if (!(MRIC->G[i])) { MRIStepCoupling_Free(MRIC); return(NULL); }
    }

    /* allocate columns of each matrix in G */
    for (i=0; i<nmat; i++)
      for (j=0; j<stages; j++) {
        MRIC->G[i][j] = NULL;
        MRIC->G[i][j] = (realtype *) calloc( stages, sizeof(realtype) );
        if (!(MRIC->G[i][j])) { MRIStepCoupling_Free(MRIC); return(NULL); }
      }
  }

  return(MRIC);
}


/*---------------------------------------------------------------
  Routine to allocate and fill a MRIStepCoupling structure
  ---------------------------------------------------------------*/
MRIStepCoupling MRIStepCoupling_Create(int nmat, int stages, int q, int p,
                                       realtype *W, realtype *G, realtype *c)
{
  int i, j, k;
  MRISTEP_METHOD_TYPE type;
  MRIStepCoupling MRIC = NULL;

  /* Check for legal inputs */
  if (nmat < 1 || stages < 1 || !c) return(NULL);

  /* Check for method coefficients and set method type */
  if (W && G)
    type = MRISTEP_IMEX;
  else if (W && !G)
    type = MRISTEP_EXPLICIT;
  else if (!W && G)
    type = MRISTEP_IMPLICIT;
  else
    return(NULL);

  /* Allocate MRIStepCoupling structure */
  MRIC = MRIStepCoupling_Alloc(nmat, stages, type);
  if (!MRIC) return(NULL);

  /* -------------------------
   * Copy the inputs into MRIC
   * ------------------------- */

  /* Method and embedding order */
  MRIC->q = q;
  MRIC->p = p;

  /* Abscissae */
  for (i=0; i<stages; i++)
    MRIC->c[i] = c[i];

  /* Coupling coefficients stored as 1D arrays of length nmat * stages * stages,
     with each stages * stages matrix stored in C (row-major) order */
  if (type == MRISTEP_EXPLICIT || type == MRISTEP_IMEX) {
    for (k = 0; k < nmat; k++)
      for (i = 0; i < stages; i++)
        for (j = 0; j < stages; j++)
          MRIC->W[k][i][j] = W[stages * (stages * k + i) + j];
  }
  if (type == MRISTEP_IMPLICIT || type == MRISTEP_IMEX) {
    for (k = 0; k < nmat; k++)
      for (i = 0; i < stages; i++)
        for (j = 0; j < stages; j++)
          MRIC->G[k][i][j] = G[stages * (stages * k + i) + j];
  }

  return(MRIC);
}


/*---------------------------------------------------------------
  Construct the MRI coupling matrix for an MIS method based on
  a given 'slow' Butcher table.
  ---------------------------------------------------------------*/
MRIStepCoupling MRIStepCoupling_MIStoMRI(ARKodeButcherTable B,
                                         int q, int p)
{
  int i, j, stages;
  booleantype padding;
  realtype Asum;
  realtype ***C;
  MRISTEP_METHOD_TYPE type;
  MRIStepCoupling MRIC;

  const realtype tol = RCONST(100.0) * UNIT_ROUNDOFF;

  /* Check that input table is non-NULL */
  if (!B) return(NULL);

  /* -----------------------------------
   * Check that the input table is valid
   * ----------------------------------- */

  /* First stage is just old solution */
  Asum = SUNRabs(B->c[0]);
  for (j=0; j<B->stages; j++)
    Asum += SUNRabs(B->A[0][j]);
  if (Asum > tol) return(NULL);

  /* Last stage exceeds 1 */
  if (B->c[B->stages-1] > ONE+tol) return(NULL);

  /* All stages are sorted */
  for (j=1; j<B->stages; j++)
    if ((B->c[j] - B->c[j-1]) < -tol) return(NULL);

  /* Each stage at most diagonally implicit */
  Asum = ZERO;
  for (i=0; i<B->stages; i++)
    for (j=i+1; j<B->stages; j++)
      Asum += SUNRabs(B->A[i][j]);
  if (Asum > tol) return(NULL);

  /* -----------------------------------------
   * determine whether the table needs padding
   * ----------------------------------------- */

  padding = SUNFALSE;

  /* Last stage time should equal 1 */
  if (SUNRabs(B->c[B->stages-1] - ONE) > tol)
    padding = SUNTRUE;

  /* Last row of A should equal b */
  for (j=0; j<B->stages; j++) {
    if (SUNRabs(B->A[B->stages-1][j] - B->b[j]) > tol)
      padding = SUNTRUE;
  }
  stages = (padding) ? B->stages+1 : B->stages;

  /* -------------------------
   * determine the method type
   * ------------------------- */

  /* Check if the table is strictly lower triangular (explicit) */
  type = MRISTEP_EXPLICIT;

  for (i=0; i<B->stages; i++)
    for (j=i; j<B->stages; j++)
      if (SUNRabs(B->A[i][j]) > tol)
        type = MRISTEP_IMPLICIT;

  /* ----------------------------
   * construct coupling structure
   * ---------------------------- */

  MRIC = MRIStepCoupling_Alloc(1, stages, type);
  if (!MRIC) return(NULL);

  /* Copy method/embedding orders */
  MRIC->q = q;
  MRIC->p = p;

  /* Copy abscissae, padding if needed */
  for (i=0; i<B->stages; i++)
    MRIC->c[i] = B->c[i];

  if (padding)
    MRIC->c[stages-1] = ONE;

  /* Construct the coupling table */
  if (type == MRISTEP_EXPLICIT)
    C = MRIC->W;
  else
    C = MRIC->G;

  /* First row is identically zero */
  for (i=0; i<stages; i++)
    for (j=0; j<stages; j++)
      C[0][i][j] = ZERO;

  /* Remaining rows = A(2:end,:) - A(1:end-1,:) */
  for (i=1; i<B->stages; i++)
    for (j=0; j<B->stages; j++)
      C[0][i][j] = B->A[i][j] - B->A[i-1][j];

  /* Padded row = b(:) - A(end,:) */
  if (padding)
    for (j=0; j<B->stages; j++)
      C[0][stages-1][j] = B->b[j] - B->A[B->stages-1][j];

  return(MRIC);
}


/*---------------------------------------------------------------
  Routine to copy a MRIStepCoupling structure
  ---------------------------------------------------------------*/
MRIStepCoupling MRIStepCoupling_Copy(MRIStepCoupling MRIC)
{
  int i, j, k, nmat, stages;
  MRISTEP_METHOD_TYPE type;
  MRIStepCoupling MRICcopy;

  /* Check for legal input */
  if (!MRIC) return(NULL);

  /* Check for method coefficients and set method type */
  if (MRIC->W && MRIC->G)
    type = MRISTEP_IMEX;
  else if (MRIC->W && !(MRIC->G))
    type = MRISTEP_EXPLICIT;
  else if (!(MRIC->W) && MRIC->G)
    type = MRISTEP_IMPLICIT;
  else
    return(NULL);

  /* Check for stage times */
  if (!(MRIC->c)) return(NULL);

  /* Get the number of coupling matrices and stages */
  nmat   = MRIC->nmat;
  stages = MRIC->stages;

  /* Allocate coupling structure */
  MRICcopy = MRIStepCoupling_Alloc(nmat, stages, type);
  if (!MRICcopy) return(NULL);

  /* Copy method and embedding orders */
  MRICcopy->q = MRIC->q;
  MRICcopy->p = MRIC->p;

  /* Copy abscissae */
  for (i=0; i<stages; i++)
    MRICcopy->c[i] = MRIC->c[i];

  /* Copy explicit coupling matrices W */
  if (MRIC->W)
    for (k = 0; k < nmat; k++)
      for (i = 0; i < stages; i++)
        for (j = 0; j < stages; j++)
          MRICcopy->W[k][i][j] = MRIC->W[k][i][j];

  /* Copy implicit coupling matrices G */
  if (MRIC->G)
    for (k = 0; k < nmat; k++)
      for (i = 0; i < stages; i++)
        for (j = 0; j < stages; j++)
          MRICcopy->G[k][i][j] = MRIC->G[k][i][j];

  return(MRICcopy);
}


/*---------------------------------------------------------------
  Routine to query the MRIStepCoupling structure workspace size
  ---------------------------------------------------------------*/
void MRIStepCoupling_Space(MRIStepCoupling MRIC, sunindextype *liw,
                           sunindextype *lrw)
{
  /* initialize outputs and return if MRIC is not allocated */
  *liw = 0;
  *lrw = 0;
  if (!MRIC) return;

  /* fill outputs based on MRIC */
  *liw = 4;
  if (MRIC->c)
    *lrw += MRIC->stages;
  if (MRIC->W)
    *lrw += MRIC->nmat * MRIC->stages * MRIC->stages;
  if (MRIC->G)
    *lrw += MRIC->nmat * MRIC->stages * MRIC->stages;
}


/*---------------------------------------------------------------
  Routine to free a MRIStepCoupling structure
  ---------------------------------------------------------------*/
void MRIStepCoupling_Free(MRIStepCoupling MRIC)
{
  int k, i;

  /* Free each field within MRIStepCoupling structure, and then
     free structure itself */
  if (MRIC) {

    if (MRIC->c)
      free(MRIC->c);

    if (MRIC->W) {
      for (k=0; k<MRIC->nmat; k++)
        if (MRIC->W[k]) {
          for (i=0; i<MRIC->stages; i++)
            if (MRIC->W[k][i]) {
              free(MRIC->W[k][i]);
              MRIC->W[k][i] = NULL;
            }
          free(MRIC->W[k]);
          MRIC->W[k] = NULL;
        }
      free(MRIC->W);
    }

    if (MRIC->G) {
      for (k=0; k<MRIC->nmat; k++)
        if (MRIC->G[k]) {
          for (i=0; i<MRIC->stages; i++)
            if (MRIC->G[k][i]) {
              free(MRIC->G[k][i]);
              MRIC->G[k][i] = NULL;
            }
          free(MRIC->G[k]);
          MRIC->G[k] = NULL;
        }
      free(MRIC->G);
    }

    free(MRIC);
  }
}


/*---------------------------------------------------------------
  Routine to print a MRIStepCoupling structure
  ---------------------------------------------------------------*/
void MRIStepCoupling_Write(MRIStepCoupling MRIC, FILE *outfile)
{
  int i, j, k;

  /* check for vaild coupling structure */
  if (!MRIC) return;
  if (!(MRIC->G)) return;

  if (MRIC->W) {
    for (i = 0; i < MRIC->nmat; i++) {
      if (!(MRIC->W[i])) return;
      for (j = 0; j < MRIC->stages; j++)
        if (!(MRIC->W[i][j])) return;
    }
  }

  if (MRIC->G) {
    for (i = 0; i < MRIC->nmat; i++) {
      if (!(MRIC->G[i])) return;
      for (j = 0; j < MRIC->stages; j++)
        if (!(MRIC->G[i][j])) return;
    }
  }

  if (!(MRIC->c)) return;

  fprintf(outfile, "  nmat = %i\n", MRIC->nmat);
  fprintf(outfile, "  stages = %i\n", MRIC->stages);
  fprintf(outfile, "  method order (q) = %i\n", MRIC->q);
  fprintf(outfile, "  embedding order (p) = %i\n", MRIC->p);

  fprintf(outfile, "  c = ");
  for (i = 0; i < MRIC->stages; i++)
    fprintf(outfile, "%"RSYM"  ", MRIC->c[i]);
  fprintf(outfile, "\n");

  if (MRIC->W) {
    for (k = 0; k < MRIC->nmat; k++) {
      fprintf(outfile, "  W[%i] = \n", k);
      for (i = 0; i < MRIC->stages; i++){
        fprintf(outfile, "      ");
        for (j = 0; j < MRIC->stages; j++)
          fprintf(outfile, "%"RSYMW"  ", MRIC->W[k][i][j]);
        fprintf(outfile, "\n");
      }
      fprintf(outfile, "\n");
    }
  }

  if (MRIC->G) {
    for (k = 0; k < MRIC->nmat; k++) {
      fprintf(outfile, "  G[%i] = \n", k);
      for (i = 0; i < MRIC->stages; i++) {
        fprintf(outfile, "      ");
        for (j = 0; j < MRIC->stages; j++)
          fprintf(outfile, "%"RSYMW"  ", MRIC->G[k][i][j]);
        fprintf(outfile, "\n");
      }
      fprintf(outfile, "\n");
    }
  }
}


/* ===========================================================================
 * Private Functions
 * ===========================================================================*/


/* ---------------------------------------------------------------------------
 * Stage type identifier: returns one of the constants
 *
 * MRISTAGE_ERK_FAST    -- standard MIS-like stage
 * MRISTAGE_ERK_NOFAST  -- standard ERK stage
 * MRISTAGE_DIRK_NOFAST -- standard DIRK stage
 * MRISTAGE_DIRK_FAST   -- coupled DIRK + MIS-like stage
 *
 * for each nontrivial stage in an MRI-like method. Otherwise (i.e., stage is
 * not in [1,MRIC->stages-1]), returns ARK_INVALID_TABLE (<0).
 *
 * The stage type is determined by 2 factors:
 * (a) Sum |MRIC->G[:][is][is]| (nonzero => DIRK)
 * (b) MRIC->c[is] - MRIC->c[is-1]  (nonzero => fast)
 * ---------------------------------------------------------------------------*/

int mriStepCoupling_GetStageType(MRIStepCoupling MRIC, int is)
{
  int i;
  realtype Gabs, cdiff;
  const realtype tol = RCONST(100.0) * UNIT_ROUNDOFF;

  if ((is < 1) || (is >= MRIC->stages)) return ARK_INVALID_TABLE;

  /* sum of stage diagonal entries across implicit tables */
  Gabs = ZERO;
  if (MRIC->G)
    for (i = 0; i < MRIC->nmat; i++)
      Gabs += SUNRabs(MRIC->G[i][is][is]);

  /* abscissae difference */
  cdiff = MRIC->c[is] - MRIC->c[is-1];

  if (Gabs > tol) {     /* DIRK */
    if (cdiff > tol) {  /* Fast */
      return(MRISTAGE_DIRK_FAST);
    } else {
      return(MRISTAGE_DIRK_NOFAST);
    }
  } else {              /* ERK */
    if (cdiff > tol) {  /* Fast */
      return(MRISTAGE_ERK_FAST);
    } else {
      return(MRISTAGE_ERK_NOFAST);
    }
  }
}


/* ---------------------------------------------------------------------------
 * Computes the stage RHS vector storage maps. With repeated abscissae the
 * first stage of the pair generally corresponds to a column of zeros and so
 * does not need to be computed and stored. The stage_map indicate if the RHS
 * needs to be computed and where to store it i.e., stage_map[i] > -1.
 * ---------------------------------------------------------------------------*/

int mriStepCoupling_GetStageMap(MRIStepCoupling MRIC,
                                int* stage_map,
                                int* nstages_stored)
{
  int i, j, k, idx;
  realtype Wsum, Gsum;
  const realtype tol = RCONST(100.0) * UNIT_ROUNDOFF;

  /* ----------------------
   * Check for valid inputs
   * ---------------------- */

  if (!MRIC) return(ARK_ILL_INPUT);
  if (!(MRIC->W) && !(MRIC->G)) return(ARK_ILL_INPUT);
  if (!stage_map || !nstages_stored) return(ARK_ILL_INPUT);

  /* -------------------
   * Compute storage map
   * ------------------- */

  /* Number of stage RHS vectors stored */
  *nstages_stored = 0;

  /* Initial storage index */
  idx = 0;

  /* Check if a stage corresponds to a column of zeros for all coupling
   * matrices by computing the column sums */
  for (j = 0; j < MRIC->stages; j++) {

    Wsum = ZERO;
    Gsum = ZERO;

    if (MRIC->W)
      for (k = 0; k < MRIC->nmat; k++)
        for (i = 0; i < MRIC->stages; i++)
          Wsum += SUNRabs(MRIC->W[k][i][j]);

    if (MRIC->G)
      for (k = 0; k < MRIC->nmat; k++)
        for (i = 0; i < MRIC->stages; i++)
          Gsum += SUNRabs(MRIC->G[k][i][j]);

    if (Wsum > tol || Gsum > tol) {
      stage_map[j] = idx;
      idx++;
    } else {
      stage_map[j] = -1;
    }
  }

  /* Check and set number of stage RHS vectors stored */
  if (idx < 1) return(ARK_ILL_INPUT);

  *nstages_stored = idx;

  return(ARK_SUCCESS);
}


/*===============================================================
  EOF
  ===============================================================*/
