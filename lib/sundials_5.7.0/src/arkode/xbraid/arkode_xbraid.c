/* --------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * --------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * --------------------------------------------------------------------------
 * This is the implementation file for the ARKODE + XBraid interface.
 * -------------------------------------------------------------------------- */

#include "arkode/arkode_xbraid.h"
#include "sundials/sundials_math.h"
#include "arkode/arkode.h"

#include "arkode_xbraid_impl.h"
#include "arkode_arkstep_impl.h"

#define ONE RCONST(1.0)


/* -------------------------------
 * Construct, initialize, and free
 * ------------------------------- */


/* Create XBraid app strucutre */
int ARKBraid_Create(void *arkode_mem, braid_App *app)
{
  int             flag;
  ARKBraidContent content;

  /* Check input */
  if (arkode_mem == NULL) return SUNBRAID_ILLINPUT;

  /* Create XBraid interface object */
  flag = SUNBraidApp_NewEmpty(app);
  if (flag != SUNBRAID_SUCCESS) return flag;

  /* Set operations */
  (*app)->ops->getvectmpl = ARKBraid_GetVecTmpl;

  /* Create ARKODE interface content */
  content = NULL;
  content = (ARKBraidContent) malloc(sizeof(*content));
  if (content == NULL)
  {
    (void) SUNBraidApp_FreeEmpty(app);
    return SUNBRAID_ALLOCFAIL;
  }

  /* Initialize content */

  /* Attach ARKODE memory */
  content->ark_mem = (ARKodeMem) arkode_mem;

  /* Interface functions */
  content->step   = ARKBraid_Step;
  content->init   = ARKBraid_Init;
  content->snorm  = SUNBraidVector_SpatialNorm;
  content->access = ARKBraid_Access;

  /* Saved return flags */
  content->last_flag_braid  = SUNBRAID_SUCCESS;
  content->last_flag_arkode = SUNBRAID_SUCCESS;

  /* Output time and solution (allocaed in access if necessary) */
  content->tout = content->ark_mem->tn;
  content->yout = NULL;

  /* Attach content */
  (*app)->content = content;

  return SUNBRAID_SUCCESS;
}


/* Initialize XBraid, attach interface functions */
int ARKBraid_BraidInit(MPI_Comm comm_w, MPI_Comm comm_t, realtype tstart,
                       realtype tstop, sunindextype ntime, braid_App app,
                       braid_Core *core)
{
  braid_Int       braid_flag;
  ARKBraidContent content;

  /* Check inputs */
  if (comm_w == MPI_COMM_NULL || comm_t == MPI_COMM_NULL || ntime < 2 ||
      app == NULL)
    return SUNBRAID_ILLINPUT;

  if (app->content == NULL) return SUNBRAID_MEMFAIL;

  /* Shortcut to content */
  content = (ARKBraidContent) app->content;

  /* Initialize XBraid */
  braid_flag = braid_Init(comm_w, comm_t, tstart, tstop, ntime, app,
                          content->step,
                          content->init,
                          SUNBraidVector_Clone,
                          SUNBraidVector_Free,
                          SUNBraidVector_Sum,
                          content->snorm,
                          content->access,
                          SUNBraidVector_BufSize,
                          SUNBraidVector_BufPack,
                          SUNBraidVector_BufUnpack,
                          core);
  CHECK_BRAID_RETURN(content->last_flag_braid, braid_flag);

  return SUNBRAID_SUCCESS;
}


/* Deallocate XBraid app structure */
int ARKBraid_Free(braid_App *app)
{
  ARKBraidContent content;  /* ARKBraid app content  */

  if (*app == NULL) return SUNBRAID_SUCCESS;

  if ((*app)->content != NULL)
  {
    content = (ARKBraidContent) (*app)->content;

    if (content->yout != NULL)
    {
      N_VDestroy(content->yout);
      content->yout = NULL;
    }
    free((*app)->content);
    (*app)->content = NULL;
  }
  return SUNBraidApp_FreeEmpty(app);
}


/* ----------------------
 * ARKBraid Set Functions
 * ---------------------- */


int ARKBraid_SetStepFn(braid_App app, braid_PtFcnStep step)
{
  ARKBraidContent content;

  if (app == NULL) return SUNBRAID_ILLINPUT;
  if (app->content == NULL) return SUNBRAID_MEMFAIL;

  content = (ARKBraidContent) app->content;

  /* Restore default or set function pointer */
  if (step == NULL)
    content->step = ARKBraid_Step;
  else
    content->step = step;

  return SUNBRAID_SUCCESS;
}


int ARKBraid_SetInitFn(braid_App app, braid_PtFcnInit init)
{
  ARKBraidContent content;

  if (app == NULL) return SUNBRAID_ILLINPUT;
  if (app->content == NULL) return SUNBRAID_MEMFAIL;

  content = (ARKBraidContent) app->content;

  /* Restore default or set function pointer */
  if (init == NULL)
    content->init = ARKBraid_Init;
  else
    content->init = init;

  return SUNBRAID_SUCCESS;
}


int ARKBraid_SetSpatialNormFn(braid_App app, braid_PtFcnSpatialNorm snorm)
{
  ARKBraidContent content;

  if (app == NULL) return SUNBRAID_ILLINPUT;
  if (app->content == NULL) return SUNBRAID_MEMFAIL;

  content = (ARKBraidContent) app->content;

  /* Restore default or set function pointer */
  if (snorm == NULL)
    content->snorm = SUNBraidVector_SpatialNorm;
  else
    content->snorm = snorm;

  return SUNBRAID_SUCCESS;
}


int ARKBraid_SetAccessFn(braid_App app, braid_PtFcnAccess access)
{
  ARKBraidContent content;

  if (app == NULL) return SUNBRAID_ILLINPUT;
  if (app->content == NULL) return SUNBRAID_MEMFAIL;

  content = (ARKBraidContent) app->content;

  /* Restore default or set function pointer */
  if (access == NULL)
    content->access = ARKBraid_Access;
  else
    content->access = access;

  return SUNBRAID_SUCCESS;
}


/* ----------------------
 * ARKBraid Get Functions
 * ---------------------- */


int ARKBraid_GetVecTmpl(braid_App app, N_Vector *tmpl)
{
  ARKBraidContent content;
  if (app == NULL) return SUNBRAID_ILLINPUT;
  if (app->content == NULL) return SUNBRAID_MEMFAIL;
  content = (ARKBraidContent) app->content;
  if (content->ark_mem == NULL) return SUNBRAID_MEMFAIL;
  *tmpl = content->ark_mem->yn;
  return SUNBRAID_SUCCESS;
}


int ARKBraid_GetARKStepMem(braid_App app, void **arkode_mem)
{
  ARKBraidContent content;
  if (app == NULL) return SUNBRAID_ILLINPUT;
  if (app->content == NULL) return SUNBRAID_MEMFAIL;
  content = (ARKBraidContent) app->content;
  if (content->ark_mem == NULL) return SUNBRAID_MEMFAIL;
  *arkode_mem = (void*) content->ark_mem;
  return SUNBRAID_SUCCESS;
}


int ARKBraid_GetUserData(braid_App app, void **user_data)
{
  ARKBraidContent content;
  if (app == NULL) return SUNBRAID_ILLINPUT;
  if (app->content == NULL) return SUNBRAID_MEMFAIL;
  content = (ARKBraidContent) app->content;
  if (content->ark_mem == NULL) return SUNBRAID_MEMFAIL;
  *user_data = content->ark_mem->user_data;
  return SUNBRAID_SUCCESS;
}


int ARKBraid_GetLastBraidFlag(braid_App app, int *last_flag)
{
  ARKBraidContent content;
  if (app == NULL) return SUNBRAID_ILLINPUT;
  if (app->content == NULL) return SUNBRAID_MEMFAIL;
  content = (ARKBraidContent) app->content;
  *last_flag = content->last_flag_braid;
  return SUNBRAID_SUCCESS;
}


int ARKBraid_GetLastARKStepFlag(braid_App app, int *last_flag)
{
  ARKBraidContent content;
  if (app == NULL) return SUNBRAID_ILLINPUT;
  if (app->content == NULL) return SUNBRAID_MEMFAIL;
  content = (ARKBraidContent) app->content;
  *last_flag = content->last_flag_arkode;
  return SUNBRAID_SUCCESS;
}


int ARKBraid_GetSolution(braid_App app, realtype *tout, N_Vector yout)
{
  ARKBraidContent content;
  if (app == NULL) return SUNBRAID_ILLINPUT;
  if (app->content == NULL) return SUNBRAID_MEMFAIL;
  content = (ARKBraidContent) app->content;
  if (content->yout == NULL) return SUNBRAID_MEMFAIL;
  *tout = content->tout;
  N_VScale(ONE, content->yout, yout);
  return SUNBRAID_SUCCESS;
}


/* --------------------------
 * XBraid Interface Functions
 * -------------------------- */


/* Take a time step */
int ARKBraid_Step(braid_App app, braid_Vector ustop, braid_Vector fstop,
                  braid_Vector u, braid_StepStatus status)
{
  braid_Int       braid_flag; /* braid function return flag  */
  int             ark_flag;   /* arkode step return flag     */
  int             flag;       /* arkode function return flag */
  int             level;      /* current level               */
  int             rfac;       /* refinement factor           */
  realtype        tstart;     /* current time                */
  realtype        tstop;      /* evolve to this time         */
  realtype        hacc;       /* accuracy based step size    */
  ARKBraidContent content;    /* ARKBraid app content        */

  /* Check input */
  if (app == NULL || status == NULL) return SUNBRAID_ILLINPUT;
  if (app->content == NULL || u->y == NULL) return SUNBRAID_MEMFAIL;

  /* Access app content */
  content = (ARKBraidContent) app->content;

  if (content->ark_mem == NULL) return SUNBRAID_MEMFAIL;

  /* Get step start and stop times */
  braid_flag = braid_StepStatusGetTstartTstop(status, &tstart, &tstop);
  CHECK_BRAID_RETURN(content->last_flag_braid, braid_flag);

  /* Propagate the solution */
  flag = ARKBraid_TakeStep((void*)(content->ark_mem), tstart, tstop, u->y,
                           &ark_flag);
  CHECK_ARKODE_RETURN(content->last_flag_arkode, flag);

  /* Refine grid (XBraid will ignore if refinement is disabled) */

  /* Get current level (XBraid only accepts refinements on level 0) */
  braid_flag = braid_StepStatusGetLevel(status, &level);
  CHECK_BRAID_RETURN(content->last_flag_braid, braid_flag);

  /* Compute refinement factor */
  if (level == 0)
  {
    /* Default to no refinement */
    rfac = 1;

    /* The step failed due to solver failure or too much error */
    if (ark_flag != 0)
    {
      /* Get the suggested step size. The rfac value is given by ETACF on a
         solver failure and limited by ETAMIN on an error test failure */
      flag = ARKStepGetCurrentStep((void*)(content->ark_mem), &hacc);
      CHECK_ARKODE_RETURN(content->last_flag_arkode, flag);

      /* Set the refinement factor */
      rfac = (int)(SUNRceil((tstop - tstart) / hacc));

      /* Limit the refinement factor */
      rfac = (rfac < 1) ? 1 : rfac;
    }

    /* set the refinement factor */
    braid_flag = braid_StepStatusSetRFactor(status, rfac);
    CHECK_BRAID_RETURN(content->last_flag_braid, braid_flag);
  }

  return SUNBRAID_SUCCESS;
}


/* Create and initialize vectors */
int ARKBraid_Init(braid_App app, realtype t, braid_Vector *u_ptr)
{
  int             flag;     /* return flag          */
  N_Vector        y;        /* output N_Vector      */
  ARKBraidContent content;  /* ARKBraid app content */

  /* Check input */
  if (app == NULL) return SUNBRAID_ILLINPUT;
  if (app->content == NULL) return SUNBRAID_MEMFAIL;

  /* Access app content */
  content = (ARKBraidContent) app->content;

  if (content->ark_mem == NULL) return SUNBRAID_MEMFAIL;

  /* Create new NVector */
  y = NULL;
  y = N_VClone(content->ark_mem->yn);
  if (y == NULL) return SUNBRAID_ALLOCFAIL;

  /* Create new XBraid vector */
  flag = SUNBraidVector_New(y, u_ptr);
  if (flag != SUNBRAID_SUCCESS) return flag;

  /* Set initial solution at all time points */
  N_VScale(ONE, content->ark_mem->yn, y);

  return SUNBRAID_SUCCESS;
}


/* User access function */
int ARKBraid_Access(braid_App app, braid_Vector u,
                    braid_AccessStatus astatus)
{
  braid_Int       braid_flag;  /* braid return flag    */
  braid_Int       done;        /* braid finished flag  */
  braid_Int       ntpoints;    /* num pts on fine grid */
  braid_Int       idx;         /* time index for u     */
  braid_Real      time;        /* time value for u     */
  ARKBraidContent content;     /* ARKBraid app content  */

  /* Check input */
  if (app == NULL || u == NULL || astatus == NULL) return SUNBRAID_ILLINPUT;
  if (app->content == NULL || u->y == NULL) return SUNBRAID_MEMFAIL;

  /* Access app content */
  content = (ARKBraidContent) app->content;

  if (content->ark_mem) return SUNBRAID_MEMFAIL;

  /* Check if XBraid is done with the current simulation */
  braid_flag = braid_AccessStatusGetDone(astatus, &done);
  CHECK_BRAID_RETURN(content->last_flag_braid, braid_flag);

  if (done)
  {
    /* Get global number of points on the fine grid */
    braid_flag = braid_AccessStatusGetNTPoints(astatus, &ntpoints);
    CHECK_BRAID_RETURN(content->last_flag_braid, braid_flag);

    /* Get the time index for the vector u */
    braid_flag = braid_AccessStatusGetTIndex(astatus, &idx);
    CHECK_BRAID_RETURN(content->last_flag_braid, braid_flag);

    /* Get the time for the vector u */
    braid_flag = braid_AccessStatusGetT(astatus, &time);
    CHECK_BRAID_RETURN(content->last_flag_braid, braid_flag);

    /* Check if this is the last time point */
    if (idx == ntpoints - 1)
    {
      /* Allocate yout if necessary */
      if (content->yout == NULL)
      {
        content->yout = N_VClone(content->ark_mem->yn);
        if (content->yout == NULL) return SUNBRAID_ALLOCFAIL;
      }

      /* Save solution for output to user */
      content->tout = time;
      N_VScale(ONE, u->y, content->yout);
    }
  }

  return SUNBRAID_SUCCESS;
}


/* -----------------
 * Utility Functions
 * ----------------- */


/* Force a single step with ARKEvolve */
int ARKBraid_TakeStep(void *arkode_mem, realtype tstart, realtype tstop,
                      N_Vector y, int *ark_flag)
{
  int      flag;      /* generic return flag      */
  int      tmp_flag;  /* evolve return flag       */
  realtype tret;      /* return time              */

  /* Check inputs */
  if (arkode_mem == NULL) return ARK_MEM_NULL;
  if (y == NULL) return ARK_ILL_INPUT;

  /* Reset ARKStep state */
  flag = ARKStepReset(arkode_mem, tstart, y);
  if (flag != ARK_SUCCESS) return flag;

  /* Set the time step size */
  flag = ARKStepSetInitStep(arkode_mem, tstop - tstart);
  if (flag != ARK_SUCCESS) return flag;

  /* Ignore temporal error test result and force step to pass */
  flag = arkSetForcePass(arkode_mem, SUNTRUE);
  if (flag != ARK_SUCCESS) return flag;

  /* Take step, check flag below */
  tmp_flag = ARKStepEvolve(arkode_mem, tstop, y, &tret, ARK_ONE_STEP);

  /* Re-enable temporal error test check */
  flag = arkSetForcePass(arkode_mem, SUNFALSE);
  if (flag != ARK_SUCCESS) return flag;

  /* Check if evolve call failed */
  if (tmp_flag < 0)
  {
    *ark_flag = STEP_FAILED;
    return ARK_SUCCESS;
  }

  /* Check if temporal error test failed */
  flag = arkGetLastKFlag(arkode_mem, &tmp_flag);
  if (flag != ARK_SUCCESS) return flag;

  if (tmp_flag > 0)
  {
    *ark_flag = STEP_ADAPT;
    return ARK_SUCCESS;
  }

  /* Step was successful and passed the error test */
  *ark_flag = STEP_SUCCESS;
  return ARK_SUCCESS;
}
