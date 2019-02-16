/* -----------------------------------------------------------------
 * Programmer(s): Carol Woodward @ LLNL
 * -----------------------------------------------------------------
 * Example problem:
 * 
 * The following is a simple example problem, with the coding
 * needed for its solution by the accelerated fixed point solver in 
 * KINSOL. 
 * The problem is from chemical kinetics, and consists of solving 
 * the first time step in a Backward Euler solution for the 
 * following three rate equations:         
 *    dy1/dt = -.04*y1 + 1.e4*y2*y3
 *    dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e2*(y2)^2
 *    dy3/dt = 3.e2*(y2)^2
 * on the interval from t = 0.0 to t = 0.1, with initial
 * conditions: y1 = 1.0, y2 = y3 = 0. The problem is stiff.
 * Run statistics (optional outputs) are printed at the end.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <kinsol/kinsol.h>           /* access to KINSOL func., consts. */
#include <nvector/nvector_serial.h>  /* access to serial N_Vector       */
#include <sundials/sundials_types.h> /* defs. of realtype, sunindextype */

/* Problem Constants */

#define NEQ   3              /* number of equations  */
#define Y10   RCONST(1.0)    /* initial y components */
#define Y20   RCONST(0.0)
#define Y30   RCONST(0.0)
#define TOL   RCONST(1.e-10) /* function tolerance */
#define DSTEP RCONST(0.1)    /* Size of the single time step used */

#define PRIORS 2

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)

/* User-defined vector accessor macro: Ith */

/* This macro is defined in order to write code which exactly matches
   the mathematical problem description given above.

   Ith(v,i) references the ith component of the vector v, where i is in
   the range [1..NEQ] and NEQ is defined above. The Ith macro is defined
   using the N_VIth macro in nvector.h. N_VIth numbers the components of
   a vector starting from 0.
*/

#define Ith(v,i)    NV_Ith_S(v,i-1)       /* Ith numbers components 1..NEQ */


/* Private functions */

static int funcRoberts(N_Vector u, N_Vector f, void *user_data);
static void PrintOutput(N_Vector u);
static void PrintFinalStats(void *kmem);
static int check_flag(void *flagvalue, const char *funcname, int opt);
static int check_ans(N_Vector u, realtype rtol, realtype atol);

/*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 */

int main()
{
  realtype fnormtol, fnorm;
  N_Vector y, scale;
  int flag;
  void *kmem;

  fnorm = 0.0;
  y = scale = NULL;
  kmem = NULL;

  /* -------------------------
   * Print problem description
   * ------------------------- */
  
  printf("Example problem from chemical kinetics solving\n"); 
  printf("the first time step in a Backward Euler solution for the\n");
  printf("following three rate equations:\n");
  printf("    dy1/dt = -.04*y1 + 1.e4*y2*y3\n");
  printf("    dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e2*(y2)^2\n");
  printf("    dy3/dt = 3.e2*(y2)^2\n");
  printf("on the interval from t = 0.0 to t = 0.1, with initial\n");
  printf("conditions: y1 = 1.0, y2 = y3 = 0.\n"); 
  printf("Solution method: Anderson accelerated fixed point iteration.\n");

  /* --------------------------------------
   * Create vectors for solution and scales
   * -------------------------------------- */

  y = N_VNew_Serial(NEQ);
  if (check_flag((void *)y, "N_VNew_Serial", 0)) return(1);

  scale = N_VNew_Serial(NEQ);
  if (check_flag((void *)scale, "N_VNew_Serial", 0)) return(1);

  /* -----------------------------------------
   * Initialize and allocate memory for KINSOL
   * ----------------------------------------- */

  kmem = KINCreate();
  if (check_flag((void *)kmem, "KINCreate", 0)) return(1);

  /* y is used as a template */

  /* Set number of prior residuals used in Anderson acceleration */
  flag = KINSetMAA(kmem, PRIORS);

  flag = KINInit(kmem, funcRoberts, y);
  if (check_flag(&flag, "KINInit", 1)) return(1);

  /* -------------------
   * Set optional inputs 
   * ------------------- */

  /* Specify stopping tolerance based on residual */

  fnormtol  = TOL; 
  flag = KINSetFuncNormTol(kmem, fnormtol);
  if (check_flag(&flag, "KINSetFuncNormTol", 1)) return(1);

  /* -------------
   * Initial guess 
   * ------------- */

  N_VConst_Serial(ZERO, y);
  Ith(y,1) = ONE;

  /* ----------------------------
   * Call KINSol to solve problem 
   * ---------------------------- */

  /* No scaling used */
  N_VConst_Serial(ONE,scale);

  /* Call main solver */
  flag = KINSol(kmem,           /* KINSol memory block */
                y,              /* initial guess on input; solution vector */
                KIN_FP,         /* global strategy choice */
                scale,          /* scaling vector, for the variable cc */
                scale);         /* scaling vector for function values fval */
  if (check_flag(&flag, "KINSol", 1)) return(1);


  /* ------------------------------------
   * Print solution and solver statistics 
   * ------------------------------------ */

  /* Get scaled norm of the system function */

  flag = KINGetFuncNorm(kmem, &fnorm);
  if (check_flag(&flag, "KINGetfuncNorm", 1)) return(1);

#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("\nComputed solution (||F|| = %Lg):\n\n",fnorm);
#else
  printf("\nComputed solution (||F|| = %g):\n\n",fnorm);
#endif  
  PrintOutput(y);

  PrintFinalStats(kmem);

  /* check the solution error */
  flag = check_ans(y, RCONST(1e-4), RCONST(1e-6));

  /* -----------
   * Free memory 
   * ----------- */
  
  N_VDestroy_Serial(y);
  N_VDestroy_Serial(scale);
  KINFree(&kmem);

  return(flag);
}

/*
 *--------------------------------------------------------------------
 * PRIVATE FUNCTIONS
 *--------------------------------------------------------------------
 */

/* 
 * System function 
 */

static int funcRoberts(N_Vector y, N_Vector g, void *user_data)
{
  realtype y1, y2, y3;
  realtype yd1, yd3;

  y1 = Ith(y,1);
  y2 = Ith(y,2);
  y3 = Ith(y,3);

  yd1 = DSTEP * ( RCONST(-0.04)*y1 + RCONST(1.0e4)*y2*y3 );
  yd3 = DSTEP * RCONST(3.0e2)*y2*y2;

  Ith(g,1) = yd1 + Y10;
  Ith(g,2) = -yd1 - yd3 + Y20;
  Ith(g,3) = yd3 + Y30;

  return(0);
}

/* 
 * Print solution at selected points
 */

static void PrintOutput(N_Vector y)
{
  realtype y1, y2, y3;

  y1 = Ith(y,1);
  y2 = Ith(y,2);
  y3 = Ith(y,3);

#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("y =%14.6Le  %14.6Le  %14.6Le\n", y1, y2, y3);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("y =%14.6e  %14.6e  %14.6e\n", y1, y2, y3);
#else
  printf("y =%14.6e  %14.6e  %14.6e\n", y1, y2, y3);
#endif

  return;
}

/* 
 * Print final statistics
 */

static void PrintFinalStats(void *kmem)
{
  long int nni, nfe;
  int flag;
  
  /* Main solver statistics */

  flag = KINGetNumNonlinSolvIters(kmem, &nni);
  check_flag(&flag, "KINGetNumNonlinSolvIters", 1);
  flag = KINGetNumFuncEvals(kmem, &nfe);
  check_flag(&flag, "KINGetNumFuncEvals", 1);

  printf("\nFinal Statistics.. \n\n");
  printf("nni      = %6ld    nfe     = %6ld \n", nni, nfe);
}

/*
 * Check function return value...
 *    opt == 0 means SUNDIALS function allocates memory so check if
 *             returned NULL pointer
 *    opt == 1 means SUNDIALS function returns a flag so check if
 *             flag >= 0
 *    opt == 2 means function allocates memory so check if returned
 *             NULL pointer 
 */

static int check_flag(void *flagvalue, const char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, 
            "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1);
  }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr,
              "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
	      funcname, *errflag);
      return(1); 
    }
  }

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr,
            "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1);
  }

  return(0);
}

/* compare the solution to a reference solution computed with a
   tolerance of 1e-14 */
static int check_ans(N_Vector u, realtype rtol, realtype atol)
{
  int      passfail=0;        /* answer pass (0) or fail (1) flag */  
  N_Vector ref;               /* reference solution vector        */
  N_Vector ewt;               /* error weight vector              */
  realtype err;               /* wrms error                       */

  /* create reference solution and error weight vectors */
  ref = N_VClone(u);
  ewt = N_VClone(u);

  /* set the reference solution data */
  NV_Ith_S(ref,0) = RCONST(9.9678538655358029e-01);
  NV_Ith_S(ref,1) = RCONST(2.9530060962800345e-03);
  NV_Ith_S(ref,2) = RCONST(2.6160735013975683e-04);

  /* compute the error weight vector */
  N_VAbs(ref, ewt);
  N_VScale(rtol, ewt, ewt);
  N_VAddConst(ewt, atol, ewt);
  if (N_VMin(ewt) <= ZERO) {
    fprintf(stderr, "\nSUNDIALS_ERROR: check_ans failed - ewt <= 0\n\n");
    return(-1);
  }
  N_VInv(ewt, ewt);   

  /* compute the solution error */
  N_VLinearSum(ONE, u, -ONE, ref, ref);
  err = N_VWrmsNorm(ref, ewt);

  /* is the solution within the tolerances? */
  passfail = (err < ONE) ? 0 : 1; 

  if (passfail) {
    fprintf(stdout, "\nSUNDIALS_WARNING: check_ans error=%g \n\n", err);
  }

  /* Free vectors */
  N_VDestroy(ref);
  N_VDestroy(ewt);

  return(passfail);
}
