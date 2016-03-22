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
 * the CVKLU solver. See fcvode.h for usage.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include "fcvode.h"
#include "cvode_impl.h"
#include <cvode/cvode_klu.h>
 
/*
 * ----------------------------------------------------------------
 * Function : FCV_KLU
 * ----------------------------------------------------------------
 */

void FCV_KLU(int *neq, int *nnz, int *ordering, int *ier)
{
  *ier = CVKLU(CV_cvodemem, *neq, *nnz);
  CVKLUSetOrdering(CV_cvodemem, *ordering);
  CV_ls = CV_LS_KLU;
}

/*
 * ----------------------------------------------------------------
 * Function : FCV_KLUReinit
 * ----------------------------------------------------------------
 */

void FCV_KLUREINIT(int *neq, int *nnz, int *reinit_type, int *ier)
{
  *ier = CVKLUReInit(CV_cvodemem, *neq, *nnz, *reinit_type);
}

