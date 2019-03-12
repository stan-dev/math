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
 * Implementation header file for the Sundials interface to 
 * the KLU linear solver.
 * -----------------------------------------------------------------
 */

#ifndef _SUNKLU_IMPL_H
#define _SUNKLU_IMPL_H

#ifndef _S_KLU_H
#define _S_KLU_H
#include "klu.h"
#endif

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * -----------------------------------------------------------------
 * Definition of KLUData
 * -----------------------------------------------------------------
 */
 
typedef struct KLUDataRec {
 
  /* Structure for KLU-specific data */
 
  klu_symbolic *s_Symbolic;
  klu_numeric  *s_Numeric;
  klu_common    s_Common;
  int           s_ordering;
  int          (*sun_klu_solve)(klu_symbolic*, klu_numeric*, int, int, double*, klu_common*);
 
} *KLUData;
 
#ifdef __cplusplus
} 
#endif 
 
#endif 
