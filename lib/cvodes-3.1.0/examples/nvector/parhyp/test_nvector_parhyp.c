/* ----------------------------------------------------------------- 
 * Programmer(s): Slaven Peles @ LLNL
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
 * This is the testing routine to check the NVECTOR Parallel module 
 * implementation. 
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <sundials/sundials_types.h>
#include <nvector/nvector_parhyp.h>
#include <sundials/sundials_math.h>
#include "test_nvector.h"

#include <mpi.h>

#if defined( SUNDIALS_HAVE_POSIX_TIMERS) && defined(_POSIX_TIMERS)
#include <time.h>
#include <unistd.h>
#endif

/* default local vector length */
#define VECLEN 5000


static int Test_N_VMake(HYPRE_ParVector W, int myid);
static int Test_N_VGetVectorID(N_Vector W, int myid);

/* ----------------------------------------------------------------------
 * Main NVector Testing Routine
 * --------------------------------------------------------------------*/
int main(int argc, char *argv[]) 
{
  int              fails = 0;                   /* counter for test failures */
  int              globfails = 0;               /* counter for test failures */
  sunindextype     local_length, global_length; /* vector lengths            */
  N_Vector         W, X, Y, Z;                  /* test vectors              */
  MPI_Comm         comm;                        /* MPI Communicator          */
  int              nprocs, myid;                /* Number of procs, proc id  */
  sunindextype     veclen;                      /* vector length             */
  int              print_timing;
  HYPRE_Int       *partitioning;                /* Vector Partitioning       */
  HYPRE_ParVector  Xhyp;                        /* Instantiate hypre parallel vector */
  int              mpierr;                      /* mpi error flag            */

  /* check input and set vector length */
  if (argc < 3) {
    SetTiming(0);
  } else {
   print_timing = atoi(argv[2]);
   SetTiming(print_timing);
  }
  
  if (argc < 2) {
    veclen = VECLEN;
  } else {
    veclen = atol(argv[1]); 
    if (veclen <= 0) {
      printf("ERROR: length of vector must be a positive integer \n");
      return(-1); 
    }
  }

  /* printf("\nRunning with vector length %ld \n \n", veclen); */
  /* Get processor number and total number of processes */
  MPI_Init(&argc, &argv);

  comm = MPI_COMM_WORLD;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &myid);

  /* set partitioning */
  local_length = veclen;
  global_length = nprocs*local_length;
  if(HYPRE_AssumedPartitionCheck()) {
    partitioning = (HYPRE_Int*) malloc(2*sizeof(HYPRE_Int));
    partitioning[0] = myid*local_length;
    partitioning[1] = (myid+1)*local_length;
  } else {
    partitioning = (HYPRE_Int*) malloc((nprocs+1)*sizeof(HYPRE_Int));
    if (veclen <= 0) {
      printf("Using global partition.\n");
      printf("I don't do this stuff. Now exiting...\n");
      return -1;
    }
  }
  /* Create template hypre vector */
  HYPRE_ParVectorCreate(comm, global_length, partitioning, &Xhyp);
  HYPRE_ParVectorInitialize(Xhyp);

  /* Create empty vector */
  W = N_VNewEmpty_ParHyp(comm, local_length, global_length);

  /* NVector Test */

  /* Hypre specific tests */
  fails += Test_N_VMake(Xhyp, myid);

  /* Create hypre vector wrapper */
  X = N_VMake_ParHyp(Xhyp);

  /* Memory allocation tests */
  fails += Test_N_VCloneEmpty(X, myid);
  fails += Test_N_VClone(X, local_length, myid);
  fails += Test_N_VCloneEmptyVectorArray(5, X, myid);
  fails += Test_N_VCloneVectorArray(5, X, local_length, myid);

  /* Create a couple of more vectors by cloning X */
  Y = N_VClone(X);
  Z = N_VClone(X);
  
  /* Skipped tests */
  /* Accessing HYPRE vector raw data is not allowed from N_Vector interface */
  /* fails += Test_N_VSetArrayPointer(W, local_length, myid); */
  /* fails += Test_N_VGetArrayPointer(X, local_length, myid); */
  
  /* N_Vector interface tests */
  fails += Test_N_VGetVectorID(W, myid);
  fails += Test_N_VConst(X, local_length, myid);
  fails += Test_N_VLinearSum(X, Y, Z, local_length, myid);
  fails += Test_N_VProd(X, Y, Z, local_length, myid);
  fails += Test_N_VDiv(X, Y, Z, local_length, myid);
  fails += Test_N_VScale(X, Z, local_length, myid);
  fails += Test_N_VAbs(X, Z, local_length, myid);
  fails += Test_N_VInv(X, Z, local_length, myid);
  fails += Test_N_VAddConst(X, Z, local_length, myid);
  fails += Test_N_VDotProd(X, Y, local_length, global_length, myid);
  fails += Test_N_VMaxNorm(X, local_length, myid);
  fails += Test_N_VWrmsNorm(X, Y, local_length, myid);
  fails += Test_N_VWrmsNormMask(X, Y, Z, local_length, global_length, myid);
  fails += Test_N_VMin(X, local_length, myid);
  fails += Test_N_VWL2Norm(X, Y, local_length, global_length, myid);
  fails += Test_N_VL1Norm(X, local_length, global_length, myid);
  fails += Test_N_VCompare(X, Z, local_length, myid);
  fails += Test_N_VInvTest(X, Z, local_length, myid);
  fails += Test_N_VConstrMask(X, Y, Z, local_length, myid);
  fails += Test_N_VMinQuotient(X, Y, local_length, myid);

  /* Free vectors */
  N_VDestroy_ParHyp(W);
  N_VDestroy_ParHyp(X);
  N_VDestroy_ParHyp(Y);
  N_VDestroy_ParHyp(Z);
  
  /* Print result */
  if (fails) {
    printf("FAIL: NVector module failed %i tests, Proc %d \n \n", fails, myid);
  } else {
     if(myid == 0) {
       printf("SUCCESS: NVector module passed all tests, Proc %d \n \n",myid);
     }
  }

  /* check if any other process failed */
  mpierr = MPI_Allreduce(&fails, &globfails, 1, MPI_INT, MPI_MAX, comm);

  /* Free hypre template vector */
  HYPRE_ParVectorDestroy(Xhyp);
  
  MPI_Finalize();

  return(globfails);
}

/* ----------------------------------------------------------------------
 * Check vector
 * 
 * Checks if all elements of vector X are set to value ans. 
 * --------------------------------------------------------------------*/
int check_ans(realtype ans, N_Vector X, sunindextype local_length)
{
  int      failure = 0;
  sunindextype i;
  hypre_ParVector *Xvec = N_VGetVector_ParHyp(X);
  realtype *Xdata = Xvec == NULL ? NULL : hypre_VectorData(hypre_ParVectorLocalVector(Xvec));
  
  /* check vector data */
  for(i=0; i < local_length; i++) {
    failure += FNEQ(Xdata[i], ans);
  }

  if (failure > ZERO)
    return(1);
  else
    return(0);
}

/* ----------------------------------------------------------------------
 * has_data
 *
 * Utility that drills down to the hypre vector local data block and 
 * checks if it is allocated. Does not verify the size of the block.
 * --------------------------------------------------------------------*/
booleantype has_data(N_Vector X)
{
  hypre_ParVector *Xvec = N_VGetVector_ParHyp(X);
  realtype *Xdata = Xvec == NULL ? NULL : hypre_VectorData(hypre_ParVectorLocalVector(Xvec));
  if (Xdata == NULL)
    return SUNFALSE;
  else
    return SUNTRUE;
}


/* ----------------------------------------------------------------------
 * set_element
 *
 * Sets single element in hypre vector by accessing its raw block. 
 * Probably not the most efficient way to set the entire vector.
 * --------------------------------------------------------------------*/
void set_element(N_Vector X, sunindextype i, realtype val)
{
  hypre_ParVector *Xvec = N_VGetVector_ParHyp(X);
  realtype *Xdata = hypre_VectorData(hypre_ParVectorLocalVector(Xvec));
  Xdata[i] = val;
}
 
/* ----------------------------------------------------------------------
 * get_element
 *
 * Reads single element from hypre vector by accessing its raw block. 
 * Probably not the most efficient way to get the vector values.
 * --------------------------------------------------------------------*/
realtype get_element(N_Vector X, sunindextype i)
{
  hypre_ParVector *Xvec = N_VGetVector_ParHyp(X);
  const realtype *Xdata = hypre_VectorData(hypre_ParVectorLocalVector(Xvec));
  return Xdata[i];
}

/* ----------------------------------------------------------------------
 * N_VGetVectorID Test
 *
 * --------------------------------------------------------------------*/
int Test_N_VGetVectorID(N_Vector v, int myid)
{
  N_Vector_ID id = N_VGetVectorID(v);
  if (id == SUNDIALS_NVEC_PARHYP) {
    if (myid == 0) {
      printf("    PASSED test -- N_VGetVectorID \n");
      /* PRINT_TIME("    N_VMake Time: %22.15e \n \n", stop_time - start_time); */
    }
    return (0);
  } else {
    printf(">>> FAILED test -- N_VGetVectorID, Proc %d \n", myid);
    printf("    Unrecognized vector type %d \n \n", id);
    return (1);    
  }
}


/* ----------------------------------------------------------------------
 * N_VMake Test
 *
 * NOTE: This routine depends on N_VConst to check vector data.
 * --------------------------------------------------------------------*/
int Test_N_VMake(HYPRE_ParVector W, int myid)
{
  int failure;
  /* double   start_time, stop_time; */
  N_Vector X;
  int local_length = hypre_ParVectorLastIndex(W) 
                     - hypre_ParVectorFirstIndex(W) + 1;

  /* clone vector */
  /* start_time = get_time(); */  
  X = N_VMake_ParHyp(W);
  /* stop_time = get_time();  */

  /* check cloned vector */
  if (X == NULL) {
    printf(">>> FAILED test -- N_VMake, Proc %d \n", myid);
    printf("    After N_VMakeEmpty, X == NULL \n \n");
    return(1);
  } 

  /* check cloned vector data */
  if (!has_data(X)) {
    printf(">>> FAILED test -- N_VMake, Proc %d \n", myid);
    printf("    Vector data == NULL \n \n");
    N_VDestroy(X);
    return(1);
  }    

  N_VConst(ONE,X);
  failure = check_ans(ONE, X, local_length);
  if (failure) {
    printf(">>> FAILED test -- N_VMake, Proc %d \n", myid);
    printf("    Failed N_VConst check \n \n");
    N_VDestroy(X);
    return(1);
  }    

  N_VDestroy(X); 

  if (myid == 0) {
    printf("    PASSED test -- N_VMake \n");
    /* PRINT_TIME("    N_VMake Time: %22.15e \n \n", stop_time - start_time); */
  }

  return(0);
}

