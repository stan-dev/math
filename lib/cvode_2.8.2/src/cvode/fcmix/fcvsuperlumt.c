/*
 * -----------------------------------------------------------------
 * $Revision: 4402 $
 * $Date: 2015-02-28 19:35:39 -0800 (Sat, 28 Feb 2015) $
 * -----------------------------------------------------------------
 * Programmer(s): Carol Woodward @ LLNL
 * -----------------------------------------------------------------
 * LLNS Copyright Start
 * Copyright (c) 2015, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Lawrence Livermore National Laboratory in part under 
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS Copyright End
 * -----------------------------------------------------------------
 * This is the implementation file for the Fortran interface to
 * the CVSuperLUMT solver. See fcvode.h for usage.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include "fcvode.h"
#include "cvode_impl.h"
#include <cvode/cvode_superlumt.h>
 
/*
 * ----------------------------------------------------------------
 * Function : FCV_SUPERLUMT
 * ----------------------------------------------------------------
 */

void FCV_SUPERLUMT(int *nthreads, int *neq, int *nnz, int *ordering, int *ier)
{
  *ier = CVSuperLUMT(CV_cvodemem, *nthreads, *neq, *nnz);
  CVSuperLUMTSetOrdering(CV_cvodemem, *ordering);
  CV_ls = CV_LS_SUPERLUMT;
}


