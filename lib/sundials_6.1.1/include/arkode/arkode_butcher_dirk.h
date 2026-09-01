/* -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the header file for ARKode's built-in DIRK Butcher tables.
 * -----------------------------------------------------------------*/

#ifndef _ARKODE_DIRK_TABLES_H
#define _ARKODE_DIRK_TABLES_H

#include <arkode/arkode_butcher.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif


/* Butcher table accessor IDs
     ERK:    0 -  99
     DIRK: 100 - 199
     MRI:  200 - 299          */

/* DEPRECATED SDIRK_2_1_2: use ARKODE_SDIRK_2_1_2 */
#define SDIRK_2_1_2             100
/* DEPRECATED BILLINGTON_3_3_2: use ARKODE_BILLINGTON_3_3_2 */
#define BILLINGTON_3_3_2        101
/* DEPRECATED TRBDF2_3_3_2: use ARKODE_TRBDF2_3_3_2 */
#define TRBDF2_3_3_2            102
/* DEPRECATED KVAERNO_4_2_3: use ARKODE_KVAERNO_4_2_3 */
#define KVAERNO_4_2_3           103
/* DEPRECATED ARK324L2SA_DIRK_4_2_3: use ARKODE_ARK324L2SA_DIRK_4_2_3 */
#define ARK324L2SA_DIRK_4_2_3   104
/* DEPRECATED CASH_5_2_4: use ARKODE_CASH_5_2_4 */
#define CASH_5_2_4              105
/* DEPRECATED CASH_5_3_4: use ARKODE_CASH_5_3_4 */
#define CASH_5_3_4              106
/* DEPRECATED SDIRK_5_3_4: use ARKODE_SDIRK_5_3_4 */
#define SDIRK_5_3_4             107
/* DEPRECATED KVAERNO_5_3_4: use ARKODE_KVAERNO_5_3_4 */
#define KVAERNO_5_3_4           108
/* DEPRECATED ARK436L2SA_DIRK_6_3_4: use ARKODE_ARK436L2SA_DIRK_6_3_4 */
#define ARK436L2SA_DIRK_6_3_4   109
/* DEPRECATED KVAERNO_7_4_5: use ARKODE_KVAERNO_7_4_5 */
#define KVAERNO_7_4_5           110
/* DEPRECATED ARK548L2SA_DIRK_8_4_5: use ARKODE_ARK548L2SA_DIRK_8_4_5 */
#define ARK548L2SA_DIRK_8_4_5   111
/* DEPRECATED ARK437L2SA_DIRK_7_3_4: use ARKODE_ARK437L2SA_DIRK_7_3_4 */
#define ARK437L2SA_DIRK_7_3_4   112
/* DEPRECATED ARK548L2SAb_DIRK_8_4_5: use ARKODE_ARK548L2SAb_DIRK_8_4_5 */
#define ARK548L2SAb_DIRK_8_4_5  113

/* Utility #defines to ensure valid input IDs for DIRK tables */

/* DEPRECATED MIN_DIRK_NUM: use ARKODE_MIN_DIRK_NUM */
#define MIN_DIRK_NUM            100

/* DEPRECATED MAX_DIRK_NUM: use ARKODE_MAX_DIRK_NUM */
#define MAX_DIRK_NUM            113

typedef enum {
  ARKODE_DIRK_NONE = -1, /* ensure enum is signed int */
  ARKODE_MIN_DIRK_NUM = 100,
  ARKODE_SDIRK_2_1_2 = ARKODE_MIN_DIRK_NUM,
  ARKODE_BILLINGTON_3_3_2,
  ARKODE_TRBDF2_3_3_2,
  ARKODE_KVAERNO_4_2_3,
  ARKODE_ARK324L2SA_DIRK_4_2_3,
  ARKODE_CASH_5_2_4,
  ARKODE_CASH_5_3_4,
  ARKODE_SDIRK_5_3_4,
  ARKODE_KVAERNO_5_3_4,
  ARKODE_ARK436L2SA_DIRK_6_3_4,
  ARKODE_KVAERNO_7_4_5,
  ARKODE_ARK548L2SA_DIRK_8_4_5,
  ARKODE_ARK437L2SA_DIRK_7_3_4,
  ARKODE_ARK548L2SAb_DIRK_8_4_5,
  ARKODE_MAX_DIRK_NUM = ARKODE_ARK548L2SAb_DIRK_8_4_5
} ARKODE_DIRKTableID;

/* Accessor routine to load built-in DIRK table */
SUNDIALS_EXPORT ARKodeButcherTable ARKodeButcherTable_LoadDIRK(ARKODE_DIRKTableID imethod);


#ifdef __cplusplus
}
#endif

#endif
