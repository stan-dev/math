/*
 * -----------------------------------------------------------------
 * $Revision$
 * $Date$
 * ----------------------------------------------------------------- 
 * Programmer(s): Carol S. Woodward @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2019, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * Implementation header file for the SUNDIALS interface to the 
 * SuperLUMT linear solver.
 * -----------------------------------------------------------------
 */

#ifndef _SUNSLUMT_IMPL_H
#define _SUNSLUMT_IMPL_H

#ifndef _SLUMT_H
#define _SLUMT_H
/* #include "pdsp_defs.h" */
#include "slu_mt_ddefs.h"
#endif

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * -----------------------------------------------------------------
 * Definition of SLUMTData
 * -----------------------------------------------------------------
 */
 
typedef struct SLUMTDataRec {
 
  /* Structure for SuperLUMT-specific data */
 
  SuperMatrix *s_A, *s_AC, *s_L, *s_U, *s_B;
  Gstat_t *Gstat;
  int *perm_r, *perm_c;
  int num_threads;
  double diag_pivot_thresh; 
  superlumt_options_t *superlumt_options;

  int s_ordering;
 
} *SLUMTData;
 
#ifdef __cplusplus
} 
#endif 
 
#endif 
