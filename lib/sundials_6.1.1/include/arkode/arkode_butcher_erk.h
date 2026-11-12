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
 * This is the header file for ARKode's built-in ERK Butcher tables.
 * -----------------------------------------------------------------*/

#ifndef _ARKODE_ERK_TABLES_H
#define _ARKODE_ERK_TABLES_H

#include <arkode/arkode_butcher.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* Butcher table accessor IDs
     ERK:    0 -  99
     DIRK: 100 - 199
     MRI:  200 - 299          */

/* DEPRECATED HEUN_EULER_2_1_2: use ARKODE_HEUN_EULER_2_1_2 */
#define HEUN_EULER_2_1_2         0
/* DEPRECATED BOGACKI_SHAMPINE_4_2_3: use ARKODE_BOGACKI_SHAMPINE_4_2_3 */
#define BOGACKI_SHAMPINE_4_2_3   1
/* DEPRECATED ARK324L2SA_ERK_4_2_3: use ARKODE_ARK324L2SA_ERK_4_2_3 */
#define ARK324L2SA_ERK_4_2_3     2
/* DEPRECATED ZONNEVELD_5_3_4: use ARKODE_ZONNEVELD_5_3_4 */
#define ZONNEVELD_5_3_4          3
/* DEPRECATED ARK436L2SA_ERK_6_3_4: use ARKODE_ARK436L2SA_ERK_6_3_4 */
#define ARK436L2SA_ERK_6_3_4     4
/* DEPRECATED SAYFY_ABURUB_6_3_4: use ARKODE_SAYFY_ABURUB_6_3_4 */
#define SAYFY_ABURUB_6_3_4       5
/* DEPRECATED CASH_KARP_6_4_5: use ARKODE_CASH_KARP_6_4_5 */
#define CASH_KARP_6_4_5          6
/* DEPRECATED FEHLBERG_6_4_5: use ARKODE_FEHLBERG_6_4_5 */
#define FEHLBERG_6_4_5           7
/* DEPRECATED DORMAND_PRINCE_7_4_5: use ARKODE_DORMAND_PRINCE_7_4_5 */
#define DORMAND_PRINCE_7_4_5     8
/* DEPRECATED ARK548L2SA_ERK_8_4_5: use ARKODE_ARK548L2SA_ERK_8_4_5 */
#define ARK548L2SA_ERK_8_4_5     9
/* DEPRECATED VERNER_8_5_6: use ARKODE_VERNER_8_5_6 */
#define VERNER_8_5_6            10
/* DEPRECATED FEHLBERG_13_7_8: use ARKODE_FEHLBERG_13_7_8 */
#define FEHLBERG_13_7_8         11
/* DEPRECATED KNOTH_WOLKE_3_3: use ARKODE_KNOTH_WOLKE_3_3 */
#define KNOTH_WOLKE_3_3         12
/* DEPRECATED ARK437L2SA_ERK_7_3_4: use ARKODE_ARK437L2SA_ERK_7_3_4 */
#define ARK437L2SA_ERK_7_3_4    13
/* DEPRECATED ARK548L2SAb_ERK_8_4_5: use ARKODE_ARK548L2SAb_ERK_8_4_5 */
#define ARK548L2SAb_ERK_8_4_5   14

/* Utility #defines to ensure valid input IDs for ERK tables */

/* DEPRECATED MIN_ERK_NUM: use ARKODE_MIN_ERK_NUM */
#define MIN_ERK_NUM              0
/* DEPRECATED MAX_ERK_NUM: use ARKODE_MAX_ERK_NUM */
#define MAX_ERK_NUM             14

typedef enum {
  ARKODE_ERK_NONE = -1, /* ensure enum is signed int */
  ARKODE_MIN_ERK_NUM = 0,
  ARKODE_HEUN_EULER_2_1_2 = ARKODE_MIN_ERK_NUM,
  ARKODE_BOGACKI_SHAMPINE_4_2_3,
  ARKODE_ARK324L2SA_ERK_4_2_3,
  ARKODE_ZONNEVELD_5_3_4,
  ARKODE_ARK436L2SA_ERK_6_3_4,
  ARKODE_SAYFY_ABURUB_6_3_4,
  ARKODE_CASH_KARP_6_4_5,
  ARKODE_FEHLBERG_6_4_5,
  ARKODE_DORMAND_PRINCE_7_4_5,
  ARKODE_ARK548L2SA_ERK_8_4_5,
  ARKODE_VERNER_8_5_6,
  ARKODE_FEHLBERG_13_7_8,
  ARKODE_KNOTH_WOLKE_3_3,
  ARKODE_ARK437L2SA_ERK_7_3_4,
  ARKODE_ARK548L2SAb_ERK_8_4_5,
  ARKODE_MAX_ERK_NUM = ARKODE_ARK548L2SAb_ERK_8_4_5
} ARKODE_ERKTableID;

/* Accessor routine to load built-in ERK table */
SUNDIALS_EXPORT ARKodeButcherTable ARKodeButcherTable_LoadERK(ARKODE_ERKTableID imethod);


#ifdef __cplusplus
}
#endif

#endif
