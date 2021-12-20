/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner, Cody J. Balos @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * This demonstration problem simulates the advection and reaction of three
 * chemical species, u, v, and w, in a three dimensional domain. The reaction
 * mechanism is a variation of the Brusselator problem from chemical kinetics.
 * This is a PDE system with 3 components, Y = [u,v,w], satisfying the
 * equations,
 *
 *    u_t = -c * dot(grad,u) + A - (w+1) * u + v * u^2
 *    v_t = -c * dot(grad,v) + w * u - v * u^2
 *    w_t = -c * dot(grad,w) + (B - w) / ep - w * u
 *
 * for t in [0,tf], X = (x,y,z) where in (x,y,z) in [0,xmax] with periodic
 * boundary conditions. The initial condition is
 *
 *    u(0,X) = k1 * A / k4 + p(X)
 *    v(0,X) = k2 * k4 * B / (k1 * k3 * A) + p(X)
 *    w(0,X) = 3.0 + p(X)
 *    p(X)   = alpha * e^( -((X - mu)^T Sigma^{-1} (x-mu)) / (2*sqrt(|Sigma|*(2pi)^3)) )
 *
 * alpha = 0.1, mu = (xmax/2.0, xmax/2.0, xmax/2.0), and Sigma = diag(xmax/4.0).
 * The reaction rates are set so k_1 = k_2 = k_3 = k_4 = k, and k_5 = k_6
 * = 1/5e-6. The spatial derivatives are discretized with first-order upwind
 * finite differences. NOUT outputs are printed at equal intervals, and run
 * statistics are printed at the end.
 *
 * Command line options:
 *   --help             prints this message
 *   --dont-save        do not save the solution to the filesystem at the nout interval (default is to save)
 *   --output-dir       the directory where all output files will be written
 *   --nout <int>       number of output times
 *   --method           ERK, ARK-DIRK, ARK-IMEX (default), CV-BDF, CV-ADAMS, IDA
 *   --nls              nonlinear solver to use; options are newton,
 *                      tl-newton (task-local newton), or fixedpoint
 *   --fpaccel          the number of fixed-point acceleration vectors to use
 *                      (only valid when using fixedpoint nonlinear solver)
 *   --nopre            turn off preconditioning
 *   --order <int>      the method order to use
 *   --npts <int>       number of mesh points in each direction
 *   --xmax <realtype>  maximum value of x (size of domain)
 *   --tf <realtype>    final time
 *   --A <realtype>     A parameter value
 *   --B <realtype>     B parameter value
 *   --k <realtype>     reaction rate
 *   --c <realtype>     advection speed
 *   --rtol <realtype>  relative tolerance
 *   --atol <realtype>  absolute tolerance
 * --------------------------------------------------------------------------*/

#include "advection_reaction_3D.hpp"

/* Main Program */
int main(int argc, char *argv[])
{
  SUNContext ctx;

  /* Initialize MPI */
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Init(&argc, &argv);

  /* Create SUNDIALS context */
  SUNContext_Create((void*) &comm, &ctx);

  /* Create SUNDIALS memory helper */
#if defined(USE_CUDA)
  SUNMemoryHelper mem_helper = SUNMemoryHelper_Cuda(ctx);
#elif defined(USE_HIP)
  SUNMemoryHelper mem_helper = SUNMemoryHelper_Hip(ctx);
#else
  SUNMemoryHelper mem_helper = SUNMemoryHelper_Sys(ctx);
#endif

  {
    /* general problem variables */
    N_Vector     y = NULL;      /* empty solution vector        */
    UserData     udata(ctx);    /* user data                    */
    UserOptions  uopt;          /* user options                 */
    int          retval;        /* reusable error-checking flag */
    char         fname[MXSTR];

    SUNDIALS_CXX_MARK_FUNCTION(udata.prof);

    /* Process input args and setup the problem */
    retval = SetupProblem(argc, argv, &udata, &uopt, mem_helper, ctx);
    if (check_retval(&retval, "SetupProblem", 1, udata.myid)) MPI_Abort(comm, 1);

    /* Create solution vector */
    y = N_VMake_MPIPlusX(udata.comm, LocalNvector(udata.grid->neq, ctx), ctx);
    if (check_retval((void *) y, "N_VMake_MPIPlusX", 0, udata.myid)) MPI_Abort(comm, 1);

    /* Enabled fused vector ops */
    if (uopt.fused)
    {
      retval = EnableFusedVectorOps(y);
      if (check_retval(&retval, "EnableFusedVectorOps", 1, udata.myid)) MPI_Abort(comm, 1);
    }

    /* Set the initial condition */
    retval = SetIC(y, &udata);
    if (check_retval(&retval, "SetIC", 1, udata.myid)) MPI_Abort(comm, 1);

    /* Output spatial mesh to disk (add extra point for periodic BC) */
    if (udata.myid == 0 && uopt.nout > 0)
    {
      snprintf(fname, MXSTR, "%s/mesh.txt", uopt.outputdir);
      udata.grid->MeshToFile(fname);
    }

    /* Integrate in time */
    if (uopt.method == "ERK")           retval = EvolveProblemExplicit(y, &udata, &uopt);
    else if (uopt.method == "ARK-DIRK") retval = EvolveProblemDIRK(y, &udata, &uopt);
    else if (uopt.method == "ARK-IMEX") retval = EvolveProblemIMEX(y, &udata, &uopt);
    else if (uopt.method == "CV-BDF")   retval = EvolveProblemBDF(y, &udata, &uopt);
    else if (uopt.method == "CV-ADAMS") retval = EvolveProblemAdams(y, &udata, &uopt);
    else if (uopt.method == "IDA")      retval = EvolveDAEProblem(y, &udata, &uopt);

    if (check_retval(&retval, "Evolve", 1, udata.myid)) MPI_Abort(comm, 1);

    /* Clean up */
    N_VDestroy(N_VGetLocalVector_MPIPlusX(y));
    N_VDestroy(y);
  }

  SUNMemoryHelper_Destroy(mem_helper);
  SUNContext_Free(&ctx);
  MPI_Finalize();
  return(0);
}


/* Destructor for problem data */
UserData::~UserData()
{
  /* free solution masks */
  N_VDestroy(N_VGetLocalVector_MPIPlusX(umask));
  N_VDestroy(umask);
  N_VDestroy(vmask);
  N_VDestroy(wmask);

  /* free the parallel grid */
  delete grid;

  /* close output streams */
  if (uopt->nout > 0)
  {
    if (UFID) fclose(UFID);
    if (VFID) fclose(VFID);
    if (WFID) fclose(WFID);
    if (TFID && myid == 0) fclose(TFID);
  }
}


/* --------------------------------------------------------------
 * Communication functions
 * --------------------------------------------------------------*/

/* Exchanges the boundary conditions only, */
int ExchangeBCOnly(N_Vector y, UserData* udata)
{
  int ierr;
  MPI_Status stat;
  MPI_Request reqR, reqS;

  /* shortcuts */
  int nvar  = udata->grid->dof;
  int myid  = udata->myid;
  int first = 0;
  int last  = udata->nprocs - 1;

  /* extract the data */
  realtype* Ydata = GetVecData(y);
  realtype* Wsend = udata->grid->getSendBuffer("WEST");

  /* open the East Irecv buffer */
  if (myid == last)
  {
    ierr = MPI_Irecv(udata->grid->getRecvBuffer("EAST"), nvar, MPI_SUNREALTYPE, first,
                     MPI_ANY_TAG, udata->comm, &reqR);
  }

  /* send first mesh node to the last processor */
  if (myid == first)
  {
    RAJA::forall< EXEC_POLICY >( RAJA::RangeSegment(0, nvar),
      [=] DEVICE_FUNC (int var) {
      Wsend[IDX(nvar, 0, var)] = Ydata[IDX(nvar, 0, var)];
    });
    ierr = MPI_Isend(Wsend, nvar, MPI_SUNREALTYPE,
                     last, 0, udata->comm, &reqS);
  }

  if (myid == last)
  {
    /* wait for exchange to finish */
    ierr = MPI_Wait(&reqR, &stat);
    if (ierr != MPI_SUCCESS)
    {
      fprintf(stderr, "\nERROR: error in MPI_Wait = %d\n", ierr);
      return -1;
    }
  }

  if (myid == first)
  {
    /* wait for exchange to finish */
    ierr = MPI_Wait(&reqS, &stat);
    if (ierr != MPI_SUCCESS)
    {
      fprintf(stderr, "\nERROR: error in MPI_Wait = %d\n", ierr);
      return -1;
    }
  }

  return(0);
}


/* Starts the exchange of the neighbor information */
int ExchangeAllStart(N_Vector y, UserData* udata)
{
  SUNDIALS_MARK_BEGIN(udata->prof, "Neighbor Exchange");

  /* shortcuts */
  realtype c = udata->c;

  /* extract the data */
  RAJA::View<realtype, RAJA::Layout<NDIMS+1> > Yview(GetVecData(y),
                                                     udata->grid->nxl,
                                                     udata->grid->nyl,
                                                     udata->grid->nzl,
                                                     udata->grid->dof);

  if (c > 0.0)
  {
    /* Flow moving in the positive directions uses backward difference. */
    udata->grid->ExchangeStart(
      [=] (realtype*, realtype* Esend, realtype*, realtype* Nsend, realtype* Bsend, realtype*) {
        int nxl = udata->grid->nxl;
        int nyl = udata->grid->nyl;
        int nzl = udata->grid->nzl;
        int dof = udata->grid->dof;

        auto range = RAJA::make_tuple(RAJA::RangeSegment(0, std::max(1,nxl-1)),
                                      RAJA::RangeSegment(0, std::max(1,nyl-1)),
                                      RAJA::RangeSegment(0, std::max(1,nzl-1)));

        RAJA::View<realtype, RAJA::Layout<3> >
          Eview(Esend, nyl, nzl, dof);
        RAJA::View<realtype, RAJA::Layout<3> >
          Nview(Nsend, nxl, nzl, dof);
        RAJA::View<realtype, RAJA::Layout<3> >
          Bview(Bsend, nxl, nyl, dof);

        RAJA::kernel<XYZ_KERNEL_POL>(range,
          [=] DEVICE_FUNC (int i, int j, int k) {

          if (nxl > 1)
          {
            Eview(j,k,0) = Yview(nxl-1,j,k,0);
            Eview(j,k,1) = Yview(nxl-1,j,k,1);
            Eview(j,k,2) = Yview(nxl-1,j,k,2);
          }

          if (nyl > 1)
          {
            Nview(i,k,0) = Yview(i,nyl-1,k,0);
            Nview(i,k,1) = Yview(i,nyl-1,k,1);
            Nview(i,k,2) = Yview(i,nyl-1,k,2);
          }

          if (nzl > 1)
          {
            Bview(i,j,0) = Yview(i,j,nzl-1,0);
            Bview(i,j,1) = Yview(i,j,nzl-1,1);
            Bview(i,j,2) = Yview(i,j,nzl-1,2);
          }

        });
      });
  }
  else if (c < 0.0)
  {
    /* Flow moving in the negative directions uses forward difference. */

    udata->grid->ExchangeStart(
      [=] (realtype* Wsend, realtype*, realtype*Ssend, realtype*, realtype*, realtype* Fsend) {
        auto range = RAJA::make_tuple(RAJA::RangeSegment(0, udata->grid->nxl-1),
                                      RAJA::RangeSegment(0, udata->grid->nyl-1),
                                      RAJA::RangeSegment(0, udata->grid->nzl-1));

        RAJA::View<realtype, RAJA::Layout<3> >
          Wview(Wsend, udata->grid->nyl, udata->grid->nzl, udata->grid->dof);
        RAJA::View<realtype, RAJA::Layout<3> >
          Sview(Ssend, udata->grid->nxl, udata->grid->nzl, udata->grid->dof);
        RAJA::View<realtype, RAJA::Layout<3> >
          Fview(Fsend, udata->grid->nxl, udata->grid->nyl, udata->grid->dof);

        RAJA::kernel<XYZ_KERNEL_POL>(range,
          [=] DEVICE_FUNC (int i, int j, int k) {
          Wview(j,k,0) = Yview(0,j,k,0);
          Wview(j,k,1) = Yview(0,j,k,1);
          Wview(j,k,2) = Yview(0,j,k,2);

          Sview(i,k,0) = Yview(i,0,k,0);
          Sview(i,k,1) = Yview(i,0,k,1);
          Sview(i,k,2) = Yview(i,0,k,2);

          Fview(i,j,0) = Yview(i,j,0,0);
          Fview(i,j,1) = Yview(i,j,0,1);
          Fview(i,j,2) = Yview(i,j,0,2);
        });
      });
  }

  SUNDIALS_MARK_END(udata->prof, "Neighbor Exchange");
  return(0);
}


/* Completes the exchange of the neighbor information */
int ExchangeAllEnd(UserData* udata)
{
  SUNDIALS_MARK_BEGIN(udata->prof, "Neighbor Exchange");
  udata->grid->ExchangeEnd();
  SUNDIALS_MARK_END(udata->prof, "Neighbor Exchange");
  return(0);
}


/* --------------------------------------------------------------
 * Problem setup
 * --------------------------------------------------------------*/

/* Parses the CLI arguments */
int ParseArgs(int argc, char *argv[], UserData* udata, UserOptions* uopt)
{
  /* check for input args */
  if (argc > 1)
  {
    /* loop over input args and get value */
    for (int i = 1; i < argc; i++)
    {
      string argvi(argv[i]);

      if (argvi.compare("--help") == 0)
      {
        InputError(argv[0]);
        return(-1);
      }
      else if (argvi.compare("--nout") == 0)
      {
        uopt->nout = atoi(argv[++i]);
      }
      else if (argvi.compare("--dont-save") == 0)
      {
        uopt->save = 0;
      }
      else if (argvi.compare("--output-dir") == 0)
      {
        if (strlen(argv[i+1]) > MXSTR)
        {
          if (udata->myid == 0)
            fprintf(stderr, "ERROR: output directory string is too long\n");
          return(-1);
        }
        uopt->outputdir = argv[++i];
      }
      else if (argvi.compare("--npts") == 0)
      {
        uopt->npts = atoi(argv[++i]);
      }
      else if (argvi.compare("--npxyz") == 0)
      {
        uopt->npxyz[0] = atoi(argv[++i]);
        uopt->npxyz[1] = atoi(argv[++i]);
        uopt->npxyz[2] = atoi(argv[++i]);
      }
      else if (argvi.compare("--xmax") == 0)
      {
        udata->xmax = strtod(argv[++i], NULL);
      }
      else if (argvi.compare("--A") == 0)
      {
        udata->A = strtod(argv[++i], NULL);
      }
      else if (argvi.compare("--B") == 0)
      {
        udata->B = strtod(argv[++i], NULL);
      }
      else if (argvi.compare("--k") == 0)
      {
        udata->k1 = strtod(argv[++i], NULL);
        udata->k2 = strtod(argv[++i], NULL);
        udata->k3 = strtod(argv[++i], NULL);
        udata->k4 = strtod(argv[++i], NULL);
      }
      else if (argvi.compare("--c") == 0)
      {
        udata->c = strtod(argv[++i], NULL);
      }
      else if (argvi.compare("--order") == 0)
      {
        uopt->order = atoi(argv[++i]);
      }
      else if (argvi.compare("--method") == 0)
      {
        uopt->method = string(argv[++i]);
        if (uopt->method != "ERK" &&
            uopt->method != "ARK-DIRK" &&
            uopt->method != "ARK-IMEX" &&
            uopt->method != "CV-BDF" &&
            uopt->method != "CV-ADAMS" &&
            uopt->method != "IDA")
        {
          fprintf(stderr, "ERROR: unknown method\n");
          InputError(argv[0]);
          return(-1);
        }
      }
      else if (argvi.compare("--fpaccel") == 0)
      {
        uopt->fpaccel = atoi(argv[++i]);
      }
      else if (argvi.compare("--nls") == 0)
      {
        uopt->nls = string(argv[++i]);
        if (uopt->nls != "newton" &&
            uopt->nls != "tl-newton" &&
            uopt->nls != "fixedpoint" &&
            uopt->nls != "none")
        {
          fprintf(stderr, "ERROR: unknown nls\n");
          InputError(argv[0]);
          return(-1);
        }
      }
      else if (argvi.compare("--nopre") == 0)
      {
        uopt->precond = 0;
      }
      else if (argvi.compare("--fused") == 0)
      {
        uopt->fused = 1;
      }
      else if (argvi.compare("--tf") == 0)
      {
        uopt->tf = strtod(argv[++i], NULL);
      }
      else if (argvi.compare("--rtol") == 0)
      {
        uopt->rtol = strtod(argv[++i], NULL);
      }
      else if (argvi.compare("--atol") == 0)
      {
        uopt->atol = strtod(argv[++i], NULL);
      }
      else
      {
        InputError(argv[0]);
        return(-1);
      }
    }
  }

  /* Explicit method uses no nonlinear solver */
  if (uopt->method == "ERK")
    uopt->nls = "none";

  /* CV Adams method only uses fixedpoint nonlinear solver */
  if (uopt->method == "CV-ADAMS")
    uopt->nls = "fixedpoint";

  return(0);
}


/* Fills the mask vector for the component so that
   u = y .* umask, v = y .* vmask, w = y .* wmask */
int ComponentMask(N_Vector mask, int component, const UserData* udata)
{
  SUNDIALS_CXX_MARK_FUNCTION(udata->prof);

  N_VConst(0.0, mask);

  RAJA::View<realtype, RAJA::Layout<NDIMS+1> > mask_view(GetVecData(mask),
                                                         udata->grid->nxl,
                                                         udata->grid->nyl,
                                                         udata->grid->nzl,
                                                         udata->grid->dof);
  auto range = RAJA::make_tuple(RAJA::RangeSegment(0, udata->grid->nxl),
                                RAJA::RangeSegment(0, udata->grid->nyl),
                                RAJA::RangeSegment(0, udata->grid->nzl));
  RAJA::kernel<XYZ_KERNEL_POL>(range,
    [=] DEVICE_FUNC (int xi, int yi, int zi) {
    mask_view(xi,yi,zi,component) = 1.0;
  });

  return 0;
}


/* Parses the CLI arguments and sets up the problem */
int SetupProblem(int argc, char *argv[], UserData* udata, UserOptions* uopt,
                 SUNMemoryHelper memhelper, SUNContext ctx)
{
  constexpr int STENCIL_WIDTH = 1;

  SUNDIALS_CXX_MARK_FUNCTION(udata->prof);

  /* Local variables */
  int retval = 0;
  char fname[MXSTR];

  /* MPI variables */
  udata->comm = MPI_COMM_WORLD;
  MPI_Comm_rank(udata->comm, &udata->myid);
  MPI_Comm_size(udata->comm, &udata->nprocs);

  /* Default problem parameters */
  udata->add_reactions = true;
  udata->xmax  = 1.0;
  udata->A     = 1.0;
  udata->B     = 3.5;
  udata->k1    = 1.0;
  udata->k2    = 1.0;
  udata->k3    = 1.0;
  udata->k4    = 1.0;
  udata->k5    = 1.0/5.0e-6;
  udata->k6    = 1.0/5.0e-6;
  udata->c     = 0.01;
  udata->uopt  = uopt;
  udata->TFID  = NULL;
  udata->UFID  = NULL;
  udata->VFID  = NULL;
  udata->WFID  = NULL;
  udata->nnlfi = 0;

  /* Set default integrator options */
  uopt->npxyz[0]  = 0;            /* number of processesors in x */
  uopt->npxyz[1]  = 0;            /* number of processesors in y */
  uopt->npxyz[2]  = 0;            /* number of processesors in z */
  uopt->npts      = 100;          /* number of mesh points in each direction */
  uopt->order     = 3;            /* method order             */
  uopt->method    = "ARK-DIRK";   /* stepper/method           */
  uopt->t0        = 0.0;          /* initial time             */
  uopt->tf        = 10.0;         /* final time               */
  uopt->rtol      = 1.0e-6;       /* relative tolerance       */
  uopt->atol      = 1.0e-9;       /* absolute tolerance       */
  uopt->nls       = "newton";     /* default to newton        */
  uopt->fpaccel   = 3;            /* default number of fixed point acceleration vectors */
  uopt->precond   = 1;            /* by default, precondition when appropriate */
  uopt->fused     = 0;            /* use fused vector ops     */
  uopt->save      = 1;            /* save solution to disk    */
  uopt->nout      = 10;           /* number of output times   */
  uopt->outputdir = (char *) "."; /* output directory         */

  /* Parse CLI args and set udata/uopt appropriately */
  retval = ParseArgs(argc, argv, udata, uopt);
  if (check_retval((void*)&retval, "ParseArgs", 1, udata->myid)) return -1;

  /* Setup the parallel decomposition */
  const sunindextype npts[] = {uopt->npts, uopt->npts, uopt->npts};
  const realtype amax[] = {0.0, 0.0, 0.0};
  const realtype bmax[] = {udata->xmax, udata->xmax, udata->xmax};
  udata->grid = new ParallelGrid<realtype,sunindextype,NDIMS>(memhelper,
    &udata->comm, amax, bmax, npts, 3, BoundaryType::PERIODIC, StencilType::UPWIND, STENCIL_WIDTH, uopt->npxyz
  );

  /* Create the solution masks */
  udata->umask = N_VMake_MPIPlusX(udata->comm, LocalNvector(udata->grid->neq, ctx), ctx);
  udata->vmask = N_VClone(udata->umask);
  udata->wmask = N_VClone(udata->umask);
  ComponentMask(udata->umask, 0, udata);
  ComponentMask(udata->vmask, 1, udata);
  ComponentMask(udata->wmask, 2, udata);

  /* Open output files for results */
  if (uopt->save)
  {
    if (udata->myid == 0)
    {
      sprintf(fname, "%s/t.%06d.txt", uopt->outputdir, udata->myid);
      udata->TFID = fopen(fname, "w");
    }

    sprintf(fname, "%s/u.%06d.txt", uopt->outputdir, udata->myid);
    udata->UFID = fopen(fname, "w");

    sprintf(fname, "%s/v.%06d.txt", uopt->outputdir, udata->myid);
    udata->VFID = fopen(fname, "w");

    sprintf(fname, "%s/w.%06d.txt", uopt->outputdir, udata->myid);
    udata->WFID = fopen(fname, "w");
  }

  /* Print problem setup */
  if (udata->myid == 0)
  {
    printf("\n\t\tAdvection-Reaction Test Problem\n\n");
    printf("Using the %s NVECTOR\n", NVECTOR_ID_STRING);
    printf("Number of Processors = %li\n", (long int) udata->nprocs);
    udata->grid->PrintInfo();
    printf("Problem Parameters:\n");
    printf("  A = %g\n", udata->A);
    printf("  B = %g\n", udata->B);
    printf("  k = %g\n", udata->k1);
    printf("  c = %g\n", udata->c);
    printf("Integrator Options:\n");
    printf("  order            = %d\n", uopt->order);
    printf("  method           = %s\n", uopt->method.c_str());
    printf("  nonlinear solver = %s\n", uopt->nls.c_str());
    printf("  fpaccel          = %d\n", uopt->fpaccel);
    printf("  preconditioner   = %d\n", uopt->precond);
    printf("  fused vector ops = %d\n", uopt->fused);
    printf("  t0               = %g\n", uopt->t0);
    printf("  tf               = %g\n", uopt->tf);
    printf("  reltol           = %.1e\n", uopt->rtol);
    printf("  abstol           = %.1e\n", uopt->atol);
    printf("  nout             = %d\n", uopt->nout);
    printf("Output directory: %s\n", uopt->outputdir);
  }


  /* return success */
  return(0);
}


/* Compute the 3D Gaussian function. */
DEVICE_FUNC
void Gaussian3D(realtype& x, realtype& y, realtype& z, realtype xmax)
{
  /* Gaussian distribution defaults */
  const realtype alpha = 0.1;
  const realtype mu[3] = { xmax/RCONST(2.0), xmax/RCONST(2.0), xmax/RCONST(2.0) };
  const realtype sigma[3] = { xmax/RCONST(4.0), xmax/RCONST(4.0), xmax/RCONST(4.0) }; // Sigma = diag(sigma)

  /* denominator = 2*sqrt(|Sigma|*(2pi)^3) */
  const realtype denom = 2.0 * sqrt((sigma[0]*sigma[1]*sigma[2])*pow(2*M_PI,3));
  x = alpha * exp( -((x - mu[0])*(x - mu[0])*(1.0/sigma[0])) / denom );
  y = alpha * exp( -((y - mu[1])*(y - mu[1])*(1.0/sigma[1])) / denom );
  z = alpha * exp( -((z - mu[2])*(z - mu[2])*(1.0/sigma[2])) / denom );
}


/* Initial condition function */
int SetIC(N_Vector y, UserData* udata)
{
  SUNDIALS_CXX_MARK_FUNCTION(udata->prof);

  /* Variable shortcuts */
  const int      nxl  = udata->grid->nxl;
  const int      nyl  = udata->grid->nyl;
  const int      nzl  = udata->grid->nzl;
  const realtype dx   = udata->grid->dx;
  const realtype dy   = udata->grid->dy;
  const realtype dz   = udata->grid->dz;
  const realtype xmax = udata->xmax;
  const realtype A    = udata->A;
  const realtype B    = udata->B;
  const realtype k1   = udata->k1;
  const realtype k2   = udata->k2;
  const realtype k3   = udata->k3;
  const realtype k4   = udata->k4;
  const int      xcrd = udata->grid->coords[0];
  const int      ycrd = udata->grid->coords[1];
  const int      zcrd = udata->grid->coords[2];

  /* Steady state solution */
  const realtype us = k1 * A / k4;
  const realtype vs = k2 * k4 * B / (k1 * k3 * A);
  const realtype ws = 3.0;

  /* Gaussian perturbation of the steady state solution */
  RAJA::View<realtype, RAJA::Layout<NDIMS+1> > yview(GetVecData(y), nxl, nyl, nzl,
                                                     udata->grid->dof);
  auto range = RAJA::make_tuple(RAJA::RangeSegment(0, nxl),
                                RAJA::RangeSegment(0, nyl),
                                RAJA::RangeSegment(0, nzl));
  RAJA::kernel<XYZ_KERNEL_POL>(range,
    [=] DEVICE_FUNC (int xi, int yi, int zi) {
    realtype x = (xcrd * nxl + xi) * dx;
    realtype y = (ycrd * nyl + yi) * dy;
    realtype z = (zcrd * nzl + zi) * dz;
    Gaussian3D(x,y,z,xmax);
    const realtype p = x + y + z;
    yview(xi,yi,zi,0) = us + p;
    yview(xi,yi,zi,1) = vs + p;
    yview(xi,yi,zi,2) = ws + p;
  });

  /* Return success */
  return(0);
}


/* Write time and solution to disk */
int WriteOutput(realtype t, N_Vector y, UserData* udata, UserOptions* uopt)
{
  SUNDIALS_CXX_MARK_FUNCTION(udata->prof);
  
  realtype  u, v, w, N;
  realtype* ydata = NULL;

  /* get vector data array */
  ydata = N_VGetArrayPointer(y);
  if (check_retval((void *) ydata, "N_VGetArrayPointer", 0, udata->myid)) return -1;

  CopyVecFromDevice(N_VGetLocalVector_MPIPlusX(y));

  /* output current solution norm to screen */
  N = (realtype) udata->grid->npts();
  u = N_VWL2Norm(y, udata->umask);
  u = sqrt(u*u/N);
  v = N_VWL2Norm(y, udata->vmask);
  v = sqrt(v*v/N);
  w = N_VWL2Norm(y, udata->wmask);
  w = sqrt(w*w/N);
  if (udata->myid == 0) {
    printf("     %10.6f   %10.6f   %10.6f   %10.6f\n", t, u, v, w);
    std::fflush(stdout);
  }

  if (uopt->save)
  {
    /* output the times to disk */
    if (udata->myid == 0 && udata->TFID)
      fprintf(udata->TFID," %.16e\n", t);

    /* output results to disk */
    RAJA::View<realtype, RAJA::Layout<NDIMS+1> > Yview(ydata,
                                                       udata->grid->nxl,
                                                       udata->grid->nyl,
                                                       udata->grid->nzl,
                                                       udata->grid->dof);

    auto range = RAJA::make_tuple(RAJA::RangeSegment(0, udata->grid->nxl),
                                  RAJA::RangeSegment(0, udata->grid->nyl),
                                  RAJA::RangeSegment(0, udata->grid->nzl));

    RAJA::kernel<XYZ_KERNEL_SERIAL_POLICY>(range,
      [=] (int i, int j, int k) {
      fprintf(udata->UFID," %.16e", Yview(i,j,k,0));
      fprintf(udata->VFID," %.16e", Yview(i,j,k,1));
      fprintf(udata->WFID," %.16e", Yview(i,j,k,2));
    });

    fprintf(udata->UFID,"\n");
    fprintf(udata->VFID,"\n");
    fprintf(udata->WFID,"\n");
  }
  
  return(0);
}


void InputError(char *name)
{
  int myid;

  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  if (myid == 0)
  {
    fprintf(stderr, "\nERROR: Invalid command line input\n");
    fprintf(stderr, "\nCommand line options for %s\n",name);
    fprintf(stderr, "  --help                    prints this message\n");
    fprintf(stderr, "  --output-dir              the directory where all output files will be written (default is the CWD)\n");
    fprintf(stderr, "  --nout <int>              number of output times to print (default is 10)\n");
    fprintf(stderr, "  --dont-save               do not save the solution to the filesystem at the nout interval (default is to save)\n");
    fprintf(stderr, "  --method                  ERK, ARK-DIRK, ARK-IMEX (default), CV-BDF, CV-ADAMS, IDA\n");
    fprintf(stderr, "  --fpaccel                 the number of fixed-point acceleration vectors to use (only valid when using fixedpoint nonlinear solver)\n");
    fprintf(stderr, "  --nls                     nonlinear solver to use (newton, tl-newton (task-local newton), fixedpoint)\n");
    fprintf(stderr, "  --nopre                   do not precondition the linear system\n");
    fprintf(stderr, "  --order <int>             the method order to use\n");
    fprintf(stderr, "  --npts <int>              number of mesh points in each direction\n");
    fprintf(stderr, "  --npxyz <int> <int> <int> number of processors in each direction (0 forces MPI to decide)\n");
    fprintf(stderr, "  --xmax <realtype>         maximum value of x (size of domain)\n");
    fprintf(stderr, "  --tf <realtype>           final time\n");
    fprintf(stderr, "  --A <realtype>            A parameter value\n");
    fprintf(stderr, "  --B <realtype>            B parameter value\n");
    fprintf(stderr, "  --k <realtype>            reaction rate\n");
    fprintf(stderr, "  --c <realtype>            advection speed\n");
    fprintf(stderr, "  --rtol <realtype>         relative tolerance\n");
    fprintf(stderr, "  --atol <realtype>         absolute tolerance\n");
  }

  MPI_Barrier(MPI_COMM_WORLD);
}

