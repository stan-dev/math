/*----------------------------------------------------------------- 
 * Programmer(s): Daniel R. Reynolds @ SMU
 *                Aaron Collier @ LLNL
 *-----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2019, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 *-----------------------------------------------------------------
 * This is the Fortran interface include file for the BBD
 * preconditioner (IDABBDPRE)
 *-----------------------------------------------------------------*/

/*==============================================================================
                    FIDABBD Interface Package
 
  The FIDABBD Interface Package is a package of C functions which, together with 
  the FIDA Interface Package, support the use of the IDA solver and MPI-parallel 
  N_Vector module, along with the IDABBDPRE preconditioner module, for the 
  solution of DAE systems in a mixed Fortran/C setting.  The combination of IDA
  and IDABBDPRE solves the linear systems arising from the solution of DAE 
  systems using a Krylov iterative linear solver via the IDASPILS interface, 
  and with a preconditioner that is block-diagonal with banded blocks.  While
  IDA and IDABBDPRE are written in C, it is assumed here that the user's
  calling program and user-supplied problem-defining routines are written in
  Fortran.
 
  The user-callable functions in this package, with the corresponding
  IDA and IDABBDPRE functions, are as follows: 

    Fortran              IDA
    --------------       ---------------------------
    FIDABBDININT         IDABBDPrecInit
    FIDABBDREINIT        IDABBDPrecReInit
    FIDABBDOPT           (accesses optional outputs)
    FIDABBDFREE          IDABBDPrecFree
    --------------       ---------------------------
 
  In addition to the Fortran residual function FIDARESFUN, the
  user-supplied functions used by this package, are listed below,
  each with the corresponding interface function which calls it (and its
  type within IDABBDPRE or IDA):

   Fortran           IDA               Type
   --------------    -----------       -----------------
   FIDAGLOCFN        FIDAgloc          IDABBDLocalFn
   FIDACOMMFN        FIDAcfn           IDABBDCommFn
   FIDAJTSETUP(*)    FIDAJTSetup       IDASpilsJTSetupFn
   FIDAJTIMES(*)     FIDAJtimes        IDASpilsJacTimesVecFn
   --------------    -----------       -----------------
   (*) = optional

  Important notes on portability:

  The names of all user-supplied routines here are fixed, in 
  order to maximize portability for the resulting mixed-language 
  program.

  Additionally, the names of the interface functions, and the names of
  the Fortran user routines called by them, appear as dummy names
  which are mapped to actual values by a series of definitions in the
  header file fidabbd.h.
 
  ==============================================================================
 
                Usage of the FIDA/FIDABBD Interface Packages
 
  The usage of the combined interface packages FIDA and FIDABBD requires
  calls to several interface functions, and a few different user-supplied
  routines which define the problem to be solved and indirectly define
  the preconditioner.  These function calls and user routines are
  summarized separately below.
 
  Some details are omitted, and the user is referred to the IDA user 
  document for more complete information.
 
  (1) User-supplied residual routine: FIDARESFUN

      The user must in all cases supply the following Fortran routine

        SUBROUTINE FIDARESFUN(T, Y, YP, R, IPAR, RPAR, IER)

      It must set the R array to F(t,y,y'), the residual of the DAE 
      system, as a function of T, Y and YP.

      The arguments are:
        T    -- scalar value of the independent variable [realtype, input]
        Y    -- array containing state variables [realtype, input]
        YP   -- array containing state derivatives [realtype, input]
        R    -- array containing DAE residuals [realtype, output]
        IPAR -- array containing integer user data that was passed
                to FIDAMALLOC [long int, input]
        RPAR -- array containing real user data that was passed to
                FIDAMALLOC [realtype, input]
        IER  -- return flag [int, output]:
                   0 if successful, 
                  >0 if a recoverable error occurred,
                  <0 if an unrecoverable error ocurred.
 
  (2) User-supplied routines to define preconditoner: FIDAGLOCFN 
      and FIDACOMMFN
 
      The routines in the IDABBDPRE module provide a preconditioner matrix
      for IDA that is block-diagonal with banded blocks.  The blocking
      corresponds to the distribution of the dependent variable vectors y 
      and y' among the processes.  Each preconditioner block is generated 
      from the Jacobian of the local part (associated with the current 
      process) of a given function G(t,y,y') approximating F(t,y,y').  The 
      blocks are generated by a difference quotient scheme independently 
      by each process, utilizing an assumed banded structure with given 
      half-bandwidths.  A separate pair of half-bandwidths defines the 
      band matrix retained.
 
  (2.1) Local approximate function FIDAGLOCFN.

      The user must supply a subroutine of the form

        SUBROUTINE FIDAGLOCFN(NLOC, T, YLOC, YPLOC, GLOC, IPAR, RPAR, IER)

      Computes the function G(t,y,y') which approximates the residual
      function F(t,y,y').  This function is to be computed locally, i.e., 
      without interprocess communication.  (The case where G is 
      mathematically identical to F is allowed.)  


      The arguments are:
        NLOC -- local problem size [long int, input]
        T    -- current time [realtype, input]
        YLOC -- array containing local state variables 
                [realtype, input]
       YPLOC -- array containing local state variable derivatives
                [realtype, input]
        GLOC -- array containing local DAE residuals [realtype, output]
        IPAR -- array containing integer user data that was passed
                to FIDAMALLOC [long int, input]
        RPAR -- array containing real user data that was passed to
                FIDAMALLOC [realtype, input]
        IER  -- return flag [int, output]:
                   0 if successful, 
                  >0 if a recoverable error occurred,
                  <0 if an unrecoverable error ocurred.
 
  (2.2) Communication function FIDACOMMF.

      The user must also supply a subroutine of the form

        SUBROUTINE FIDACOMMFN(NLOC, T, YLOC, YPLOC, IPAR, RPAR, IER)

      Performs all interprocess communication necessary to evaluate the 
      approximate residual function G described above.  It is expected to
      save communicated data in  work space defined by the user, and made 
      available to FIDAGLOCFN.  Each call to the FIDACOMMFN is preceded 
      by a call to FIDARESFUN with the same (t,y,y') arguments.  Thus 
      FIDACOMMFN can omit any communications done by FIDARESFUN if 
      relevant to the evaluation of G.

      The arguments are:
        NLOC -- local problem size [long int, input]
        T    -- current time [realtype, input]
        YLOC -- array containing local state variables 
                [realtype, input]
        YPLOC -- array containing local state variable derivatives
                [realtype, input]
        IPAR -- array containing integer user data that was passed
                to FIDAMALLOC [long int, input]
        RPAR -- array containing real user data that was passed to
                FIDAMALLOC [realtype, input]
        IER  -- return flag [int, output]:
                   0 if successful, 
                  >0 if a recoverable error occurred,
                  <0 if an unrecoverable error ocurred.

 
  (3) Optional user-supplied Jacobian-vector setup and product 
      functions: FIDAJTSETUP and FIDAJTIMSE

      As an option, the user may supply a routine that computes the 
      product of the system Jacobian J = dF/dy and a given vector v.  
      If supplied, a 'setup' routine to prepare any user data 
      structures must exist, and have the form:
 
        SUBROUTINE FIDAJTSETUP(T, Y, YP, R, CJ, EWT, H, IPAR, RPAR, IER)

      It must perform any relevant preparations for subsequent calls to 
      the user-provided FIDAJTIMES routine (see below).  

      The arguments are:
        T    -- current time [realtype, input]
        Y    -- array containing state variables [realtype, input]
        YP   -- array containing state variable derivatives
                [realtype, input]
        R    -- array containing DAE residuals [realtype, input]
        CJ   -- current value of scalar in Jacobian [realtype, input]
        EWT  -- array containing error weight vector [realtype, input]
        H    -- current step size [realtype, input]
        IPAR -- array containing integer user data that was passed to
                FIDAMALLOC [long int, input]
        RPAR -- array containing real user data that was passed to
                FIDAMALLOC [realtype, input]
        IER  -- return flag [int, output]:
                   0 if successful, 
                   nonzero if an error.
 
      The accompanying Jacobian matrix-vector product routine must 
      have the following form:

        SUBROUTINE FIDAJTIMES(T, Y, YP, R, V, FJV, CJ, EWT, H, 
       1                      IPAR, RPAR, WK1, WK2, IER)

      This routine must compute the product vector J*v, and store 
      the product in FJV.  
 
      The arguments are:
        T    -- current time [realtype, input]
        Y    -- array containing state variables [realtype, input]
        YP   -- array containing state variable derivatives
                [realtype, input]
        R    -- array containing DAE residuals [realtype, input]
        V    -- vector to multiply [realtype, input]
        FJV  -- product vector [realtype, output]
        CJ   -- current value of scalar in Jacobian [realtype, input]
        EWT  -- array containing error weight vector [realtype, input]
        H    -- current step size [realtype, input]
        IPAR -- array containing integer user data that was passed to
                FIDAMALLOC [long int, input]
        RPAR -- array containing real user data that was passed to
                FIDAMALLOC [realtype, input]
        WK1, WK2 -- arrays containing temporary workspace of same size
                as Y [realtype, input]
        IER  -- return flag [int, output]:
                    0 if successful, 
                    nonzero if an error.
 
  (4) Initialization:  FNVINITP, generic iterative linear solver initialization, 
      FIDAMALLOC, FIDASPILSINIT, and FIDABBDINIT.
 
  (4.1) To initialize the parallel machine environment, the user must make 
      the following call:

         CALL FNVINITP(COMM, 2, NLOCAL, NGLOBAL, IER)

      where the second argument is an int containing the IDA
      solver ID (2). The other arguments are:
        COMM = the MPI communicator [int, input]
        NLOCAL = local vector size on this processor 
           [long int, input]
        NGLOBAL = system size, and the global size of vectors 
           (the sum of all values of NLOCAL) [long int, input]
        IER = return completion flag [int, ouptut]. 
                  0 = success, 
                 -1 = failure.

      NOTE: The COMM argument passed to the FNVINITP routine is only supported 
      if the MPI implementation used to build SUNDIALS includes the MPI_Comm_f2c
      function from the MPI-2 specification.  To check if the function is 
      supported look for the line "#define SUNDIALS_MPI_COMM_F2C 1" in the 
      sundials_config.h header file.
 
 (4.2) To initialize a generic iterative linear solver structure for 
      solving linear systems within the Newton solver, the user must make one 
      of the following calls:

          CALL FSUNPCGINIT(2, PRETYPE, MAXL, IER)
          CALL FSUNSPBCGSINIT(2, PRETYPE, MAXL, IER)
          CALL FSUNSPFGMRINIT(2, PRETYPE, MAXL, IER)
          CALL FSUNSPGMRINIT(2, PRETYPE, MAXL, IER)
          CALL FSUNSPTFQMRINIT(2, PRETYPE, MAXL, IER)

      In each of these, one argument is an int containing the IDA solver 
      ID (2). 

      The other arguments are:

        PRETYPE = type of preconditioning to perform (0=none, 1=left, 
           2=right, 3=both) [int, input]
        MAXL = maximum Krylov subspace dimension [int, input]
        IER = return completion flag [int, output]:
	          0 = success, 
		 -1 = failure.

  (4.3) To set various problem and solution parameters and allocate
      internal memory, make the following call:

        CALL FIDAMALLOC(T0, Y0, YP0, IATOL, RTOL, ATOL, ID, CONSTR,
       1                IOUT, ROUT, IPAR, RPAR, IER)

      The arguments are:
        T0    = initial value of t [realtype, input]
        Y0    = array of initial conditions, y(t0) [realtype, input]
        YP0   = value of y'(t0) [realtype, input]
        IATOL = type for absolute tolerance ATOL [int, input]: 
                  1 = scalar, 
                  2 = array,
                  3 = user-supplied function; the user must 
                      supply a routine FIDAEWT to compute the 
		      error weight vector.
        RTOL  = scalar relative tolerance [realtype, input]
        ATOL  = scalar/array absolute tolerance [realtype, input]
        IOUT  = array of length at least 21 for integer optional 
                inputs and outputs [long int, output]
        ROUT  = array of length 6 for real optional inputs and 
                outputs [realtype, output]
	IPAR = array of user integer data [long int, in/out]
	RPAR = array with user real data [realtype, in/out]
	IER  = return completion flag [int, output]:
                  0 = SUCCESS,
                 -1 = failure (see printed message for details).
 
      The optional outputs are:
 
            LENRW   = IOUT( 1) -> IDAGetWorkSpace
            LENIW   = IOUT( 2) -> IDAGetWorkSpace
            NST     = IOUT( 3) -> IDAGetNumSteps
            NRE     = IOUT( 4) -> IDAGetNumResEvals
            NETF    = IOUT( 5) -> IDAGetNumErrTestFails
            NCFN    = IOUT( 6) -> IDAGetNumNonlinSolvConvFails
            NNI     = IOUT( 7) -> IDAGetNumNonlinSolvIters
            NSETUPS = IOUT( 8) -> IDAGetNumLinSolvSetups
            KLAST   = IOUT( 9) -> IDAGetLastOrder
            KCUR    = IOUT(10) -> IDAGetCurrentOrder
            NBCKTRK = IOUT(11) -> IDAGetNumBacktrackOps
            NGE     = IOUT(12) -> IDAGetNumGEvals
 
            HINUSED = ROUT( 1) -> IDAGetActualInitStep
            HLAST   = ROUT( 2) -> IDAGetLastStep
            HCUR    = ROUT( 3) -> IDAGetCurrentStep
            TCUR    = ROUT( 4) -> IDAGetCurrentTime
            TOLSFAC = ROUT( 5) -> IDAGetTolScaleFactor
            UNITRND = ROUT( 6) -> UNIT_ROUNDOFF

      The user data arrays IPAR and RPAR are passed unmodified to 
      all subsequent calls to user-provided routines. Changes to
      either array inside a user-provided routine will be 
      propagated. Using these two arrays, the user can dispense 
      with COMMON blocks to pass data betwen user-provided 
      routines. 

      If the user program includes the FIDAEWT routine for the 
      evaluation of the error weights, the following call must be made

         CALL FIDAEWTSET (FLAG, IER)

      with FLAG = 1 to specify that FIDAEWT is provided.
      The return flag IER is 0 if successful, and nonzero otherwise.
 
  (4.4) Create the IDASPILS interface to attach the generic 
     iterative linear solver to IDA, by making the following call:
    
       CALL FIDASPILSINIT(IER)

     The arguments are:
	IER = error return flag [int, output]: 
	       0 = success; 
	      <0 = an error occured
 
  (4.5) To allocate memory and initialize data associated with the 
      IDABBDPRE preconditioner, make the following call:

        CALL FIDABBDINIT(NLOCAL, MUDQ, MLDQ, MU, ML, DQRELY, IER)

      The arguments are:
        NLOCAL = local vector size on this process 
             [long int, input]
        MUDQ = upper half-bandwidth to be used in the computation
             of the local Jacobian blocks by difference 
             quotients.  These may be smaller than the true 
             half-bandwidths of the Jacobian of the local block 
             of g, when smaller values may provide greater 
             efficiency [long int, input]
        MLDQ = lower half-bandwidth to be used in the computation
             of the local Jacobian blocks by difference 
             quotients [long int, input]
        MU = upper half-bandwidth of the band matrix that is
             retained as an approximation of the local Jacobian
             block (may be smaller than MUDQ) [long int, input]
        ML = lower half-bandwidth of the band matrix that is
             retained as an approximation of the local Jacobian
             block (may be smaller than MLDQ) [long int, input]
        DQRELY = relative increment factor in y for difference 
             quotients [realtype, input]
                    0.0 = default (sqrt(unit roundoff))
        IER = return completion flag [int, output]:
                    0 = success
                   <0 = an error occurred
 
  (4.6) To specify whether the Krylov linear solver should use the 
      supplied FIDAJTSETUP and FIDAJTIMES routines, or the internal 
      finite difference approximation, make the call

         CALL FIDASPILSSETJAC(FLAG, IER)

      with the int FLAG=1 to specify that FIDAJTSETUP and FIDAJTIMES 
      are provided (FLAG=0 specifies to use and internal finite 
      difference approximation to this product).  The int return 
      flag IER=0 if successful, and nonzero otherwise.
 
  (5) Re-initialization: FIDAREINIT, FIDABBDREINIT

      If a sequence of problems of the same size is being solved using 
      the Krylov linear solver in combination with the IDABBDPRE 
      preconditioner, then the IDA package can be reinitialized for 
      the second and subsequent problems so as to avoid further memory 
      allocation.  First, in place of the call to FIDAMALLOC, make the 
      following call:

        CALL FIDAREINIT(T0, Y0, YP0, IATOL, RTOL, ATOL, ID, CONSTR, IER)

      The arguments have the same names and meanings as those of 
      FIDAMALLOC.  FIDAREINIT performs the same initializations as 
      FIDAMALLOC, but does no memory allocation for IDA data structures, 
      using instead the existing internal memory created by the previous 
      FIDAMALLOC call.  

      Following the call to FIDAREINIT, if there is no change in any of 
      the linear solver arguments, but the user wishes to modify the 
      values of MUDQ, MLDQ or DQRELY from the previous call to 
      FIDABBDINIT, then a user may call:

       CALL FIDABBDREINIT(MUDQ, MLDQ, DQRELY, IER)

      The arguments have the same names and meanings as those of 
      FIDABBDINIT.

      However, if there is a change in any of the linear solver 
      arguments or other preconditioner arguments, then a call to
      FSUNPCGINIT, FSUNSPBCGSINIT, FSUNSPFGMRINIT, FSUNSPGMRINIT, 
      or FSUNSPTFQMRINIT is required; in this case the linear 
      solver memory is reallocated.  Following this call, the 
      IDASPILS interface must also be reconstructed using another
      call to FIDASPILSINIT (interface memory is freed and 
      reallocated), as well as a subsequent call to FIDABBDINIT.
 
  (6) The solver: FIDASOLVE

      To solve the DAE system, make the following call:

        CALL FIDASOLVE(TOUT, TRET, Y, YP, ITASK, IER)

      The arguments are:
       TOUT = next value of t at which a solution is desired 
           [realtype, input]
       TRET = value of t reached by the solver [realtype, output]
       Y = state variable array [realtype, output]
       YP = state variable derivative array [realtype, output]
       ITASK = task indicator [int, input]:
                 1 = normal mode (overshoot TOUT and interpolate)
                 2 = one-step mode (return after each internal 
                     step taken)
                 3 = normal tstop mode (like 1, but integration 
                     never proceeds past TSTOP, which must be 
                     specified through a call to FIDASETRIN using
                     the key 'STOP_TIME')
                 4 = one step tstop (like 2, but integration 
                     never goes past TSTOP)
       IER = completion flag [int, output]: 
                  0 = success, 
                  1 = tstop return, 
                  2 = root return, 
                  values -1 ... -10 are failure modes (see 
                    IDA manual).
     The current values of the optional outputs are immediately 
     available in the IOUT and ROUT arrays.
 
  (7) Optional outputs: FIDABBDOPT

      Optional outputs specific to the IDASpils linear solver are 
      available in IOUT(13)...IOUT(21)
 
      To obtain the optional outputs associated with the IDABBDPRE 
      module, make the following call:

        CALL FIDABBDOPT (LENRWBBD, LENIWBBD, NGEBBD)

      The arguments returned are:
        LENRWBBD = length of real preconditioner work space, in 
             realtype words (this size is local to the current 
             process if run in parallel) [long int, output]
        LENIWBBD = length of integer preconditioner work space, in 
             integer words (this size is local to the current 
             process if run in parallel) [long int, output]
        NGEBBD   = number of G(t,y,y') evaluations (calls to 
             FIDAGLOCFN) so far [long int, output]
 
  (8) Memory freeing: FIDAFREE

      To the free the internal memory created by the calls to 
      FNVINITP, FIDAMALLOC, FIDASPILSINIT and FIDABBDINIT, make 
      the following call:

        CALL FIDAFREE
 
==============================================================================*/

#ifndef _FIDABBD_H
#define _FIDABBD_H

#include <sundials/sundials_nvector.h>
#include <sundials/sundials_types.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#if defined(SUNDIALS_F77_FUNC)

#define FIDA_BBDINIT    SUNDIALS_F77_FUNC(fidabbdinit, FIDABBDINIT)
#define FIDA_BBDREINIT  SUNDIALS_F77_FUNC(fidabbdreinit, FIDABBDREINIT)
#define FIDA_BBDOPT     SUNDIALS_F77_FUNC(fidabbdopt, FIDABBDOPT)
#define FIDA_GLOCFN     SUNDIALS_F77_FUNC(fidaglocfn, FIDAGLOCFN)
#define FIDA_COMMFN     SUNDIALS_F77_FUNC(fidacommfn, FIDACOMMFN)

#else

#define FIDA_BBDINIT    fidabbdinit_
#define FIDA_BBDREINIT  fidabbdreinit_
#define FIDA_BBDOPT     fidabbdopt_
#define FIDA_GLOCFN     fidaglocfn_
#define FIDA_COMMFN     fidacommfn_

#endif

/* Prototypes of exported functions */

void FIDA_BBDINIT(long int *Nloc, long int *mudq, long int *mldq,
                  long int *mu, long int *ml, realtype *dqrely, int *ier);

void FIDA_BBDREINIT(long int *Nloc, long int *mudq, long int *mldq,
		    realtype *dqrely, int *ier);

void FIDA_BBDOPT(long int *lenrwbbd, long int *leniwbbd,
                 long int *ngebbd);

/* Prototypes: Functions Called by the IDABBD Module */

int FIDAgloc(long int Nloc, realtype t, N_Vector yy, N_Vector yp,
             N_Vector gval, void *user_data);
int FIDAcfn(long int Nloc, realtype t, N_Vector yy, N_Vector yp,
            void *user_data);

#ifdef __cplusplus
}
#endif

#endif
