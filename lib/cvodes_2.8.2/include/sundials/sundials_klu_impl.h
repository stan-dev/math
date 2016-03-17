/*
 * -----------------------------------------------------------------
 * $Revision: 4075 $
 * $Date: 2014-04-24 10:46:58 -0700 (Thu, 24 Apr 2014) $
 * ----------------------------------------------------------------- 
 * Programmer(s): Carol S. Woodward @ LLNL
 * -----------------------------------------------------------------
 * LLNS Copyright Start
 * Copyright (c) 2014, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Lawrence Livermore National Laboratory in part under 
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS Copyright End
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
 
} *KLUData;
 
#ifdef __cplusplus
} 
#endif 
 
#endif 
