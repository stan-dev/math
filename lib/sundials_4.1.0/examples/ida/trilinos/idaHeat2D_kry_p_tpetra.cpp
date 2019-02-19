/* -----------------------------------------------------------------
 * Programmer(s): Slaven Peles @ LLNL
 * -----------------------------------------------------------------
 * Acknowledgements: This example is based on idaHeat2D_kry_p
 *                   example by Daniel R. Reynolds @ SMU and
 *                   Allan Taylor, Alan Hindmarsh and
 *                   Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Example problem for IDA: 2D heat equation, parallel, GMRES.
 *
 * This example solves a discretized 2D heat equation problem.
 * This version uses the Krylov solver SUNSPGMR.
 *
 * The DAE system solved is a spatial discretization of the PDE
 *          du/dt = d^2u/dx^2 + d^2u/dy^2
 * on the unit square. The boundary condition is u = 0 on all edges.
 * Initial conditions are given by u = 16 x (1 - x) y (1 - y).
 * The PDE is treated with central differences on a uniform MX x MY
 * grid. The values of u at the interior points satisfy ODEs, and
 * equations u = 0 at the boundaries are appended, to form a DAE
 * system of size N = MX * MY. Here MX = MY = 10.
 *
 * The system is actually implemented on submeshes, processor by
 * processor, with an MXSUB by MYSUB mesh on each of NPEX * NPEY
 * processors.
 *
 * The system is solved with IDA using the Krylov linear solver
 * SUNSPGMR. The preconditioner uses the diagonal elements of the
 * Jacobian only. Routines for preconditioning, required by
 * SUNSPGMR, are supplied here. The constraints u >= 0 are posed
 * for all components. Local error testing on the boundary values
 * is suppressed. Output is taken at t = 0, .01, .02, .04,
 * ..., 10.24.
 * -----------------------------------------------------------------
 * When running this example for unit testing, set the environemnt
 * variables OMP_PROC_BIND=false and OMP_NUM_THREADS=2 to run
 * without thread binding and with two OpenMP threads.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <vector>

#include <Tpetra_Core.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_Version.hpp>

#include <ida/ida.h>
#include <nvector/trilinos/SundialsTpetraVectorInterface.hpp>
#include <nvector/nvector_trilinos.h>
#include <sunlinsol/sunlinsol_spgmr.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_mpi_types.h>

#define ZERO  RCONST(0.0)
#define ONE   RCONST(1.0)
#define TWO   RCONST(2.0)

#define NOUT         11    /* Number of output times */

#define NPEX         2     /* No. PEs in x direction of PE array */
#define NPEY         2     /* No. PEs in y direction of PE array */
                           /* Total no. PEs = NPEX*NPEY */
#define MXSUB        5     /* No. x points per subgrid */
#define MYSUB        5     /* No. y points per subgrid */

/* Global spatial mesh is MX x MY = (NPEX x MXSUB) x (NPEY x MYSUB) */

/* type definitions */
typedef Sundials::TpetraVectorInterface::vector_type vector_type;
typedef vector_type::scalar_type scalar_type;
typedef vector_type::mag_type mag_type;
typedef vector_type::global_ordinal_type global_ordinal_type;
typedef vector_type::local_ordinal_type local_ordinal_type;
typedef vector_type::node_type::memory_space memory_space;
typedef vector_type::node_type::execution_space execution_space;
typedef vector_type::map_type map_type;

/* User data structure */
struct UserData {
  /* user data destructor */
  ~UserData()
  {
    if (pp != NULL)
      N_VDestroy(pp);
  }

  Teuchos::RCP<const Teuchos::Comm<int> > comm;
  int thispe, npex, npey, ixsub, jysub;
  sunindextype mx, my, mxsub, mysub;
  realtype dx, dy, coeffx, coeffy, coeffxy;
  Kokkos::View<scalar_type*, memory_space> uext; /* device array */

  Kokkos::View<scalar_type*, memory_space> send_buff_top;
  Kokkos::View<scalar_type*, memory_space> send_buff_bottom;
  Kokkos::View<scalar_type*, memory_space> send_buff_left;
  Kokkos::View<scalar_type*, memory_space> send_buff_right;
  Kokkos::View<scalar_type*, memory_space> recv_buff_top;
  Kokkos::View<scalar_type*, memory_space> recv_buff_bottom;
  Kokkos::View<scalar_type*, memory_space> recv_buff_left;
  Kokkos::View<scalar_type*, memory_space> recv_buff_right;

  typename Kokkos::View<scalar_type*, memory_space>::HostMirror h_send_buff_top;
  typename Kokkos::View<scalar_type*, memory_space>::HostMirror h_send_buff_bottom;
  typename Kokkos::View<scalar_type*, memory_space>::HostMirror h_send_buff_left;
  typename Kokkos::View<scalar_type*, memory_space>::HostMirror h_send_buff_right;
  typename Kokkos::View<scalar_type*, memory_space>::HostMirror h_recv_buff_top;
  typename Kokkos::View<scalar_type*, memory_space>::HostMirror h_recv_buff_bottom;
  typename Kokkos::View<scalar_type*, memory_space>::HostMirror h_recv_buff_left;
  typename Kokkos::View<scalar_type*, memory_space>::HostMirror h_recv_buff_right;

  N_Vector pp;    /* vector of diagonal preconditioner elements */
};


/* User-supplied residual and supporting functions */

int resHeat(realtype tt, N_Vector uu, N_Vector up,
            N_Vector rr, void *user_data);

static int rescomm(N_Vector uu, N_Vector up, void *user_data);

static int reslocal(realtype tt, N_Vector uu, N_Vector up,
                    N_Vector res,  void *user_data);

static int BSend(N_Vector uu, void *user_data);

static int BRecvPost(MPI_Request request[], void *user_data);

static int BRecvWait(MPI_Request request[], void *user_data);


/* User-supplied preconditioner functions */

int PsolveHeat(realtype tt, N_Vector uu, N_Vector up, N_Vector rr,
               N_Vector rvec, N_Vector zvec, realtype c_j,
               realtype delta, void *user_data);

int PsetupHeat(realtype tt, N_Vector yy, N_Vector yp, N_Vector rr,
               realtype c_j, void *user_data);

/* Private function to allocate memory, initialize problem and print output */

static int InitUserData(int thispe, const Teuchos::RCP<const Teuchos::Comm<int> > comm, UserData *data);

static int AllocUserData(int thispe, N_Vector uu, UserData *data);

static int SetInitialProfile(N_Vector uu, N_Vector up, N_Vector id,
                             N_Vector res, UserData *data);

static void PrintHeader(realtype rtol, realtype atol, UserData *data);

static void PrintOutput(int id, void *ida_mem, realtype t, N_Vector uu);

static void PrintFinalStats(void *ida_mem);

static int check_retval(void *flagvalue, const char *funcname, int opt, int id);

/*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 */

int main(int argc, char *argv[])
{
  void *ida_mem = NULL;
  SUNLinearSolver LS = NULL;
  int iout, retval;
  realtype rtol, atol, t0, t1, tout, tret;
  N_Vector uu = NULL;
  N_Vector up = NULL;
  N_Vector constraints = NULL;
  N_Vector id = NULL;
  N_Vector res = NULL;

  /* Start an MPI session */
  Tpetra::ScopeGuard tpetraScope(&argc, &argv);

  /* Create Tpetra communicator */
  auto comm = Tpetra::getDefaultComm();
  const int thispe = comm->getRank();
  const int npes   = comm->getSize();

  /* Allocate and initialize the data structure */
  UserData *data = new UserData();
  if(check_retval((void *)data, "malloc", 2, thispe))
    return -1;

  /* Initialize user data */
  InitUserData(thispe, comm, data);

  /* Set global solution vector lengths. */
  const sunindextype global_length = data->mx * data->my;

  /* Choose zero-based (C-style) indexing. */
  const sunindextype index_base = 0;

  /* Construct am MPI Map */
  Teuchos::RCP<const map_type> testMap =
    Teuchos::rcp(new map_type (global_length, index_base, comm,
                               Tpetra::GloballyDistributed));

  /* Check if the number of MPI processes matches the number of subgrids */
  if (npes != (data->npex * data->npey)) {
    if (thispe == 0)
      fprintf(stderr,
              "\nMPI_ERROR(0): npes = %d is not equal to NPEX*NPEY = %d\n",
              npes, data->npex * data->npey);
    delete data;
    return 1;
  }

  /* Construct a Tpetra vector and return refernce counting pointer to it. */
  Teuchos::RCP<vector_type> rcpuu =
    Teuchos::rcp(new vector_type(testMap));



  /* Allocate and initialize N-vectors. */

  uu = N_VMake_Trilinos(rcpuu);
  if(check_retval((void *)uu, "N_VMake_Trilinos", 0, thispe))
    return -1;

  up = N_VClone(uu);
  if(check_retval((void *)up, "N_VClone", 0, thispe))
    return -1;

  res = N_VClone(uu);
  if(check_retval((void *)res, "N_VClone", 0, thispe))
    return -1;

  constraints = N_VClone(uu);
  if(check_retval((void *)constraints, "N_VClone", 0, thispe))
    return -1;

  id = N_VClone(uu);
  if(check_retval((void *)id, "N_VClone", 0, thispe))
    return -1;

  /* Allocate user data extended vector and MPI buffers */
  retval = AllocUserData(thispe, uu, data);
  if(check_retval(&retval, "AllocUserData", 1, thispe)) return -1;


  /* Initialize the uu, up, id, and res profiles. */
  SetInitialProfile(uu, up, id, res, data);

  /* Set constraints to all 1's for nonnegative solution values. */
  N_VConst(ONE, constraints);

  t0 = ZERO; t1 = RCONST(0.01);

  /* Scalar relative and absolute tolerance. */
  rtol = ZERO;
  atol = RCONST(1.0e-3);

  /* Call IDACreate and IDAMalloc to initialize solution. */

  ida_mem = IDACreate();
  if(check_retval((void *)ida_mem, "IDACreate", 0, thispe))
    return -1;

  retval = IDASetUserData(ida_mem, data);
  if(check_retval(&retval, "IDASetUserData", 1, thispe))
    return -1;

  retval = IDASetSuppressAlg(ida_mem, SUNTRUE);
  if(check_retval(&retval, "IDASetSuppressAlg", 1, thispe))
    return -1;

  retval = IDASetId(ida_mem, id);
  if(check_retval(&retval, "IDASetId", 1, thispe))
    return -1;

  retval = IDASetConstraints(ida_mem, constraints);
  if(check_retval(&retval, "IDASetConstraints", 1, thispe))
    return -1;
  N_VDestroy(constraints);

  retval = IDAInit(ida_mem, resHeat, t0, uu, up);
  if(check_retval(&retval, "IDAInit", 1, thispe))
    return -1;

  retval = IDASStolerances(ida_mem, rtol, atol);
  if(check_retval(&retval, "IDASStolerances", 1, thispe))
    return -1;

  /* Call SUNLinSol_SPGMR and IDASetLinearSolver to specify the linear solver. */

  LS = SUNLinSol_SPGMR(uu, PREC_LEFT, 0);  /* use default maxl */
  if(check_retval((void *)LS, "SUNLinSol_SPGMR", 0, thispe))
    return -1;

  retval = IDASetLinearSolver(ida_mem, LS, NULL);
  if(check_retval(&retval, "IDASetLinearSolver", 1, thispe))
    return -1;

  retval = IDASetPreconditioner(ida_mem, PsetupHeat, PsolveHeat);
  if(check_retval(&retval, "IDASetPreconditioner", 1, thispe))
    return -1;

  /* Print output heading (on processor 0 only) and intial solution  */

  if (thispe == 0) PrintHeader(rtol, atol, data);
  PrintOutput(thispe, ida_mem, t0, uu);

  /* Loop over tout, call IDASolve, print output. */

  for (tout = t1, iout = 1; iout <= NOUT; iout++, tout *= TWO) {

    retval = IDASolve(ida_mem, tout, &tret, uu, up, IDA_NORMAL);
    if(check_retval(&retval, "IDASolve", 1, thispe))
      return -1;

    PrintOutput(thispe, ida_mem, tret, uu);

  }

  /* Print remaining counters. */

  if (thispe == 0) PrintFinalStats(ida_mem);

  /* Free memory */

  IDAFree(&ida_mem);
  SUNLinSolFree(LS);

  N_VDestroy(id);
  N_VDestroy(res);
  N_VDestroy(up);
  N_VDestroy(uu);

  /* Delete template Tpetra vector */
  rcpuu = Teuchos::null;

  delete data;

  return(0);
}

/*
 *--------------------------------------------------------------------
 * FUNCTIONS CALLED BY IDA
 *--------------------------------------------------------------------
 */

/*
 * resHeat: heat equation system residual function
 * This uses 5-point central differencing on the interior points, and
 * includes algebraic equations for the boundary values.
 * So for each interior point, the residual component has the form
 *    res_i = u'_i - (central difference)_i
 * while for each boundary point, it is res_i = u_i.
 *
 * This parallel implementation uses several supporting routines.
 * First a call is made to rescomm to do communication of subgrid boundary
 * data into array uext.  Then reslocal is called to compute the residual
 * on individual processors and their corresponding domains.  The routines
 * BSend, BRecvPost, and BREcvWait handle interprocessor communication
 * of uu required to calculate the residual.
 */
int resHeat(realtype tt, N_Vector uu, N_Vector up, N_Vector rr,
            void *user_data)
{
  int retval = 0;

  /* Call rescomm to do inter-processor communication. */
  retval = rescomm(uu, up, user_data);

  /* Call reslocal to calculate res. */
  retval = reslocal(tt, uu, up, rr, user_data);

  return(retval);
}


/*
 * PsetupHeat: setup for diagonal preconditioner for heatsk.
 *
 * The optional user-supplied functions PsetupHeat and
 * PsolveHeat together must define the left preconditoner
 * matrix P approximating the system Jacobian matrix
 *                   J = dF/du + cj*dF/du'
 * (where the DAE system is F(t,u,u') = 0), and solve the linear
 * systems P z = r.   This is done in this case by keeping only
 * the diagonal elements of the J matrix above, storing them as
 * inverses in a vector pp, when computed in PsetupHeat, for
 * subsequent use in PsolveHeat.
 *
 * In this instance, only cj and data (user data structure, with
 * pp etc.) are used from the PsetupHeat argument list.
 *
 */
int PsetupHeat(realtype tt, N_Vector yy, N_Vector yp, N_Vector rr,
               realtype c_j, void *user_data)
{
  /* Unwrap the user data */
  UserData* data = reinterpret_cast<UserData*>(user_data);
  const int ixsub = data->ixsub;
  const int jysub = data->jysub;
  const int npex  = data->npex;
  const int npey  = data->npey;
  const sunindextype mxsub = data->mxsub;
  const sunindextype mysub = data->mysub;

  sunindextype ibc, i0, jbc, j0;

  /* Initially set all pp elements on the device to one. */
  N_VConst(ONE, data->pp);

  /* Get access to local data */
  Teuchos::RCP<vector_type> pptp = N_VGetVector_Trilinos(data->pp);
  auto pp_2d = pptp->getLocalView<memory_space>();
  auto pp_1d = Kokkos::subview (pp_2d, Kokkos::ALL(), 0);

  /* Calculate the value for the inverse element of the diagonal preconditioner */
  const realtype pelinv = ONE/(c_j + data->coeffxy);

  ibc = (ixsub == 0) || (ixsub == npex-1);
  i0  = (ixsub == 0);
  jbc = (jysub == 0) || (jysub == npey-1);
  j0  = (jysub == 0);

  /* Set inverse of the preconditioner; ppv must be on the device */
  Kokkos::parallel_for (Kokkos::RangePolicy<execution_space>(0, (mxsub - ibc)*(mysub - jbc)),
    KOKKOS_LAMBDA (const local_ordinal_type &tid)
    {
      sunindextype j = tid / (mxsub - ibc) + j0;
      sunindextype i = tid % (mxsub - ibc) + i0;

      pp_1d(i + j*mxsub) = pelinv;
    }
  );

  return(0);
}


/*
 * PsolveHeat: solve preconditioner linear system.
 * This routine multiplies the input vector rvec by the vector pp
 * containing the inverse diagonal Jacobian elements (previously
 * computed in PsetupHeat), returning the result in zvec.
 */
int PsolveHeat(realtype tt, N_Vector uu, N_Vector up,
               N_Vector rr, N_Vector rvec, N_Vector zvec,
               realtype c_j, realtype delta, void *user_data)
{
  UserData* data = reinterpret_cast<UserData*>(user_data);

  N_VProd(data->pp, rvec, zvec);

  return(0);
}

/*
 *--------------------------------------------------------------------
 * SUPPORTING FUNCTIONS
 *--------------------------------------------------------------------
 */


/*
 * rescomm routine.  This routine performs all inter-processor
 * communication of data in u needed to calculate G.
 */
static int rescomm(N_Vector uu, N_Vector up, void* user_data)
{
  /* Set array of raw MPI requests */
  MPI_Request request[4];

  /* Start receiving boundary data from neighboring PEs. */
  BRecvPost(request, user_data);

  /* Send data from boundary of local grid to neighboring PEs. */
  BSend(uu, user_data);

  /* Finish receiving boundary data from neighboring PEs. */
  BRecvWait(request, user_data);

  return(0);
}

/*
 * reslocal routine.  Compute res = F(t, uu, up).  This routine assumes
 * that all inter-processor communication of data needed to calculate F
 * has already been done, and that this data is in the work array uext.
 */
static int reslocal(realtype tt, N_Vector uu, N_Vector up, N_Vector rr,
                    void *user_data)
{
  UserData* data = reinterpret_cast<UserData*>(user_data);

  /* Get subgrid indices, array sizes, and grid coefficients. */
  const int ixsub = data->ixsub;
  const int jysub = data->jysub;
  const int npex  = data->npex;
  const int npey  = data->npey;
  const sunindextype mxsub  = data->mxsub;
  const sunindextype mxsub2 = data->mxsub + 2;
  const sunindextype mysub  = data->mysub;
  const realtype coeffx  = data->coeffx;
  const realtype coeffy  = data->coeffy;
  const realtype coeffxy = data->coeffxy;

  sunindextype ibc, i0, jbc, j0;

  /* Initialize all elements of rr to uu. This sets the boundary
     elements simply without indexing hassles. */
  N_VScale(ONE, uu, rr);

  /* Get vector data arrays, extended work array uext. */
  Teuchos::RCP<vector_type> uuv = N_VGetVector_Trilinos(uu);
  Teuchos::RCP<vector_type> upv = N_VGetVector_Trilinos(up);
  Teuchos::RCP<vector_type> rrv = N_VGetVector_Trilinos(rr);

  const auto uu_2d = uuv->getLocalView<memory_space>();
  const auto uu_1d = Kokkos::subview (uu_2d, Kokkos::ALL(), 0);
  const auto up_2d = upv->getLocalView<memory_space>();
  const auto up_1d = Kokkos::subview (up_2d, Kokkos::ALL(), 0);
  auto rr_2d = rrv->getLocalView<memory_space>();
  auto rr_1d = Kokkos::subview (rr_2d, Kokkos::ALL(), 0);
  Kokkos::View<realtype*, memory_space> uext_1d = data->uext;

  /* Copy local segment of u vector into the working extended array uext.
     This completes uext prior to the computation of the rr vector.
     uext and uuv must be on the device.     */
  Kokkos::parallel_for (Kokkos::RangePolicy<execution_space>(0, mxsub*mysub),
    KOKKOS_LAMBDA (const local_ordinal_type &tid)
    {
      sunindextype j = tid/mxsub;
      sunindextype i = tid%mxsub;

      uext_1d((i+1) + (j+1)*mxsub2) = uu_1d(i + j*mxsub);
    }
  );

  /* Set loop limits for the interior of the local subgrid. */

  /* Prepare to loop over subgrid. */
  ibc = (ixsub == 0) || (ixsub == npex-1);
  i0  = (ixsub == 0);
  jbc = (jysub == 0) || (jysub == npey-1);
  j0  = (jysub == 0);

  /* Compute local residual; uext, upv, and resv must be on the device */
  Kokkos::parallel_for (Kokkos::RangePolicy<execution_space>(0, (mxsub - ibc)*(mysub - jbc)),
    KOKKOS_LAMBDA (const local_ordinal_type &tid)
    {
      sunindextype j = tid/(mxsub - ibc) + j0;
      sunindextype i = tid%(mxsub - ibc) + i0;
      sunindextype locu  = i + j*mxsub;
      sunindextype locue = (i+1) + (j+1)*mxsub2;

      realtype termx   = coeffx * (uext_1d(locue-1)      + uext_1d(locue+1));
      realtype termy   = coeffy * (uext_1d(locue-mxsub2) + uext_1d(locue+mxsub2));
      realtype termctr = coeffxy * uext_1d(locue);
      rr_1d(locu) = up_1d(locu) - (termx + termy - termctr);
    }
  );

  return(0);
}


/*
 * Routine to send boundary data to neighboring PEs.
 */
static int BSend(N_Vector uu, void *user_data)
{
  UserData* data = reinterpret_cast<UserData*>(user_data);

  /* Get comm, thispe, subgrid indices, data sizes */
  auto comm = Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int> >(data->comm);
  auto rawComm = comm->getRawMpiComm();

  const int thispe = data->thispe;
  const int ixsub  = data->ixsub;
  const int jysub  = data->jysub;
  const int npex   = data->npex;
  const int npey   = data->npey;

  /* Type conversion from sunindextype to int (required by Teuchos::send) */
  const int mxsub = data->mxsub;
  const int mysub = data->mysub;

  /* Get pointers to buffers and extended solution vector data array uext. */
  Kokkos::View<realtype*, memory_space> bufleft   = data->send_buff_left;
  Kokkos::View<realtype*, memory_space> bufright  = data->send_buff_right;
  Kokkos::View<realtype*, memory_space> buftop    = data->send_buff_top;
  Kokkos::View<realtype*, memory_space> bufbottom = data->send_buff_bottom;
  typename Kokkos::View<realtype*>::HostMirror h_bufleft   =  data->h_send_buff_left;
  typename Kokkos::View<realtype*>::HostMirror h_bufright  =  data->h_send_buff_right;
  typename Kokkos::View<realtype*>::HostMirror h_buftop    =  data->h_send_buff_top;
  typename Kokkos::View<realtype*>::HostMirror h_bufbottom =  data->h_send_buff_bottom;

  /* Get solution vector data. */
  Teuchos::RCP<vector_type> uuv = N_VGetVector_Trilinos(uu);
  const auto uu_2d = uuv->getLocalView<memory_space>();
  const auto uu_1d = Kokkos::subview (uu_2d, Kokkos::ALL(), 0);

  /* If jysub > 0, send data from bottom x-line of u.  (via bufbottom) */

  if (jysub != 0) {
    /* Device kernel here to copy from uarray to the buffer on the device */
    Kokkos::parallel_for (Kokkos::RangePolicy<execution_space>(0, mxsub),
      KOKKOS_LAMBDA (const local_ordinal_type &lx)
      {
        bufbottom(lx) = uu_1d(lx);
      }
    );
    /* Copy device buffer to the corresponding host buffer */
    Kokkos::deep_copy(h_bufbottom, bufbottom);
    /* MPI send buffer */
    MPI_Send(h_bufbottom.data(), mxsub, PVEC_REAL_MPI_TYPE, thispe-npex, 0, *rawComm);
  }

  /* If jysub < NPEY-1, send data from top x-line of u. (via buftop) */

  if (jysub != npey-1) {
    /* Device kernel here to copy from uarray to the buffer on the device */
    Kokkos::parallel_for (Kokkos::RangePolicy<execution_space>(0, mxsub),
      KOKKOS_LAMBDA (const local_ordinal_type &lx)
      {
        buftop(lx) = uu_1d((mysub-1)*mxsub + lx);
      }
    );
    /* Copy buffer to the host */
    Kokkos::deep_copy(h_buftop, buftop);
    /* MPI send buffer */
    MPI_Send(h_buftop.data(), mxsub, PVEC_REAL_MPI_TYPE, thispe+npex, 0, *rawComm);
  }

  /* If ixsub > 0, send data from left y-line of u (via bufleft). */

  if (ixsub != 0) {
    /* Device kernel here to copy from uarray to the buffer on the device */
    Kokkos::parallel_for (Kokkos::RangePolicy<execution_space>(0, mysub),
      KOKKOS_LAMBDA (const local_ordinal_type &ly)
      {
        bufleft(ly) = uu_1d(ly*mxsub);
      }
    );
    /* Copy buffer to the host */
    Kokkos::deep_copy(h_bufleft, bufleft);
    /* MPI send buffer */
    MPI_Send(h_bufleft.data(), mysub, PVEC_REAL_MPI_TYPE, thispe-1, 0, *rawComm);
  }

  /* If ixsub < NPEX-1, send data from right y-line of u (via bufright). */

  if (ixsub != npex-1) {
    /* Device kernel here to copy from uarray to the buffer on the device */
    Kokkos::parallel_for (Kokkos::RangePolicy<execution_space>(0, mysub),
      KOKKOS_LAMBDA (const local_ordinal_type &ly)
      {
        bufright(ly) = uu_1d(ly*mxsub + (mxsub-1));
      }
    );
    /* Copy buffer to the host */
    Kokkos::deep_copy(h_bufright, bufright);
    /* MPI send buffer */
    MPI_Send(h_bufright.data(), mysub, PVEC_REAL_MPI_TYPE, thispe+1, 0, *rawComm);
  }

  return(0);

}


/*
 * Routine to start receiving boundary data from neighboring PEs.
 * Notes:
 *   1) buffer should be able to hold 2*(MYSUB+MYSUB) realtype entries, should
 *      be passed to both the BRecvPost and BRecvWait functions, and should not
 *      be manipulated between the two calls.
 *   2) request should have 4 entries, and should be passed in
 *      both calls also.
 */
static int BRecvPost(MPI_Request request[], void *user_data)
{
  UserData* data = reinterpret_cast<UserData*>(user_data);

  /* Get comm, thispe, subgrid indices, data sizes */
  auto comm = Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int> >(data->comm);
  auto rawComm = comm->getRawMpiComm();

  const int thispe = data->thispe;
  const int ixsub  = data->ixsub;
  const int jysub  = data->jysub;
  const int npex   = data->npex;
  const int npey   = data->npey;
  const sunindextype mxsub = data->mxsub;
  const sunindextype mysub = data->mysub;

  /* Get left, right, top and bottom host buffers. */
  typename Kokkos::View<realtype*>::HostMirror h_bufleft   = data->h_recv_buff_left;
  typename Kokkos::View<realtype*>::HostMirror h_bufright  = data->h_recv_buff_right;
  typename Kokkos::View<realtype*>::HostMirror h_buftop    = data->h_recv_buff_top;
  typename Kokkos::View<realtype*>::HostMirror h_bufbottom = data->h_recv_buff_bottom;

  /* If jysub > 0, receive data for bottom x-line of uext. */
  if (jysub != 0) {
    MPI_Irecv(h_bufbottom.data(), mxsub, PVEC_REAL_MPI_TYPE,
              thispe-npex, 0, *rawComm, &request[0]);
  }

  /* If jysub < NPEY-1, receive data for top x-line of uext. */
  if (jysub != npey-1) {
    MPI_Irecv(h_buftop.data(), mxsub, PVEC_REAL_MPI_TYPE,
              thispe+npex, 0, *rawComm, &request[1]);
  }

  /* If ixsub > 0, receive data for left y-line of uext (via bufleft). */
  if (ixsub != 0) {
    MPI_Irecv(h_bufleft.data(), mysub, PVEC_REAL_MPI_TYPE,
              thispe-1, 0, *rawComm, &request[2]);
  }

  /* If ixsub < NPEX-1, receive data for right y-line of uext (via bufright). */
  if (ixsub != npex-1) {
    MPI_Irecv(h_bufright.data(), mysub, PVEC_REAL_MPI_TYPE,
              thispe+1, 0, *rawComm, &request[3]);
  }

  return(0);

}

/*
 * Routine to finish receiving boundary data from neighboring PEs.
 * Notes:
 *   1) buffer should be able to hold 2*MYSUB realtype entries, should be
 *      passed to both the BRecvPost and BRecvWait functions, and should not
 *      be manipulated between the two calls.
 *   2) request should have four entries, and should be passed in both
 *      calls also.
 */
static int BRecvWait(MPI_Request request[], void *user_data)
{
  MPI_Status status;
  UserData* data = reinterpret_cast<UserData*>(user_data);

  /* Get thispe, subgrid indices, data sizes */
  const int ixsub  = data->ixsub;
  const int jysub  = data->jysub;
  const int npex   = data->npex;
  const int npey   = data->npey;
  const sunindextype mxsub = data->mxsub;
  const sunindextype mysub = data->mysub;

  /* Get pointers to buffers and extended solution vector data array uext. */
  typename Kokkos::View<realtype*>::HostMirror h_bufleft   = data->h_recv_buff_left;
  typename Kokkos::View<realtype*>::HostMirror h_bufright  = data->h_recv_buff_right;
  typename Kokkos::View<realtype*>::HostMirror h_buftop    = data->h_recv_buff_top;
  typename Kokkos::View<realtype*>::HostMirror h_bufbottom = data->h_recv_buff_bottom;

  Kokkos::View<realtype*, memory_space> bufleft   = data->recv_buff_left;
  Kokkos::View<realtype*, memory_space> bufright  = data->recv_buff_right;
  Kokkos::View<realtype*, memory_space> buftop    = data->recv_buff_top;
  Kokkos::View<realtype*, memory_space> bufbottom = data->recv_buff_bottom;

  Kokkos::View<realtype*, memory_space> uext_1d = data->uext;

  const sunindextype mxsub2 = mxsub + 2;
  const sunindextype mysub1 = mysub + 1;

  /* If jysub > 0, receive data for bottom x-line of uext. */
  if (jysub != 0) {
    MPI_Wait(&request[0], &status);
    /* Copy the buffer from the host to the device */
    Kokkos::deep_copy(bufbottom, h_bufbottom);
    /* Copy the bottom dev_recv_buff to uext. */
    Kokkos::parallel_for (Kokkos::RangePolicy<execution_space>(0, mxsub),
      KOKKOS_LAMBDA (const local_ordinal_type &lx)
      {
        uext_1d(1 + lx) = bufbottom(lx);
      }
    );
  }

  /* If jysub < NPEY-1, receive data for top x-line of uext. */
  if (jysub != npey-1) {
    MPI_Wait(&request[1], &status);
    /* Copy the buffer from the host to the device */
    Kokkos::deep_copy(buftop, h_buftop);
    /* Copy the top dev_recv_buff to uext. */
    Kokkos::parallel_for (Kokkos::RangePolicy<execution_space>(0, mxsub),
      KOKKOS_LAMBDA (const local_ordinal_type &lx)
      {
        uext_1d((1 + mysub1*mxsub2) + lx) = buftop(lx);
      }
    );
  }

  /* If ixsub > 0, receive data for left y-line of uext (via bufleft). */
  if (ixsub != 0) {
    MPI_Wait(&request[2], &status);
    /* Copy the buffer from the host to the device */
    Kokkos::deep_copy(bufleft, h_bufleft);
    /* Copy the left dev_recv_buff to uext. */
    Kokkos::parallel_for (Kokkos::RangePolicy<execution_space>(0, mysub),
      KOKKOS_LAMBDA (const local_ordinal_type &ly)
      {
        uext_1d((ly + 1)*mxsub2) = bufleft(ly);
      }
    );
  }

  /* If ixsub < NPEX-1, receive data for right y-line of uext (via bufright). */
  if (ixsub != npex-1) {
    MPI_Wait(&request[3], &status);
    /* Copy the buffer from the host to the device */
    Kokkos::deep_copy(bufright, h_bufright);
    /* Copy the right dev_recv_buff to uext. */
    Kokkos::parallel_for (Kokkos::RangePolicy<execution_space>(0, mysub),
      KOKKOS_LAMBDA (const local_ordinal_type &ly)
      {
        uext_1d((ly + 2)*mxsub2 - 1) = bufright(ly);
      }
    );
  }

  return(0);
}

/*
 *--------------------------------------------------------------------
 * PRIVATE FUNCTIONS
 *--------------------------------------------------------------------
 */

/*
 * InitUserData initializes the user's data block data.
 */
static int InitUserData(int thispe, const Teuchos::RCP<const Teuchos::Comm<int> > comm, UserData *data)
{
  data->comm    = comm;
  data->thispe  = thispe;
  data->npex    = NPEX;  /* Number of subgrids in x-direction */
  data->npey    = NPEY;  /* Number of subgrids in y-direction */
  data->mxsub   = MXSUB; /* Number of subgrid mesh points in x-direction */
  data->mysub   = MYSUB; /* Number of subgrid mesh points in y-direction */
  data->jysub   = thispe/data->npex;
  data->ixsub   = thispe - (data->jysub * data->npex);
  data->mx      = data->npex * data->mxsub;  /* Mesh size in x-direction */
  data->my      = data->npey * data->mysub;  /* Mesh size in y-direction */
  data->dx      = ONE/(data->mx-ONE); /* Assumes a [0,1] interval in x. */
  data->dy      = ONE/(data->my-ONE); /* Assumes a [0,1] interval in y. */
  data->coeffx  = ONE/(data->dx * data->dx);
  data->coeffy  = ONE/(data->dy * data->dy);
  data->coeffxy = TWO/(data->dx * data->dx) + TWO/(data->dy * data->dy);

  return(0);
}


/*
 * AllocUserData allocates memory for the extended vector uext
 * and MPI communication buffers.
 */
static int AllocUserData(int thispe, N_Vector uu, UserData *data)
{
  sunindextype mxsub = data->mxsub;
  sunindextype mysub = data->mysub;

  /* Allocate local extended vector (includes ghost nodes) */
  Kokkos::resize(data->uext, (mxsub + 2)*(mysub +2));

  /* Allocate local send buffers */
  Kokkos::resize(data->send_buff_left,   mysub);
  Kokkos::resize(data->send_buff_right,  mysub);
  Kokkos::resize(data->send_buff_top,    mxsub);
  Kokkos::resize(data->send_buff_bottom, mxsub);

  data->h_send_buff_left = create_mirror_view(data->send_buff_left);
  data->h_send_buff_right = create_mirror_view(data->send_buff_right);
  data->h_send_buff_top = create_mirror_view(data->send_buff_top);
  data->h_send_buff_bottom = create_mirror_view(data->send_buff_bottom);

  /* Allocate local receive buffers */
  Kokkos::resize(data->recv_buff_left,   mysub);
  Kokkos::resize(data->recv_buff_right,  mysub);
  Kokkos::resize(data->recv_buff_top,    mxsub);
  Kokkos::resize(data->recv_buff_bottom, mxsub);

  data->h_recv_buff_left = create_mirror_view(data->recv_buff_left);
  data->h_recv_buff_right = create_mirror_view(data->recv_buff_right);
  data->h_recv_buff_top = create_mirror_view(data->recv_buff_top);
  data->h_recv_buff_bottom = create_mirror_view(data->recv_buff_bottom);

  /* An N-vector to hold preconditioner. */
  data->pp = N_VClone(uu);
  if(data->pp == NULL) {
    return -1;
  }

  return 0;
}


/*
 * SetInitialProfile sets the initial values for the problem.
 */
static int SetInitialProfile(N_Vector uu, N_Vector up,  N_Vector id,
                             N_Vector res, UserData *data)
{
  sunindextype i, iloc, j, jloc, loc;
  realtype xfact, yfact;

  /* Initialize uu. */

  /* Get Tpetra vectors */
  Teuchos::RCP<vector_type> tu = N_VGetVector_Trilinos(uu);
  Teuchos::RCP<vector_type> ti = N_VGetVector_Trilinos(id);

  /* Get access to local uu and id data, sync the host with the device */
  tu->sync<Kokkos::HostSpace>();
  auto u_2d = tu->getLocalView<Kokkos::HostSpace>();
  auto u_1d = Kokkos::subview(u_2d, Kokkos::ALL(), 0);
  tu->modify<Kokkos::HostSpace>();

  ti->sync<Kokkos::HostSpace>();
  auto id_2d = ti->getLocalView<Kokkos::HostSpace>();
  auto id_1d = Kokkos::subview(id_2d, Kokkos::ALL(), 0);
  ti->modify<Kokkos::HostSpace>();

  /* Set mesh spacings and subgrid indices for this PE. */
  const realtype dx = data->dx;
  const realtype dy = data->dy;
  const int ixsub = data->ixsub;
  const int jysub = data->jysub;

  /* Set beginning and ending locations in the global array corresponding
     to the portion of that array assigned to this processor. */
  const sunindextype mxsub   = data->mxsub;
  const sunindextype mysub   = data->mysub;
  const sunindextype ixbegin = mxsub*ixsub;
  const sunindextype ixend   = mxsub*(ixsub+1) - 1;
  const sunindextype jybegin = mysub*jysub;
  const sunindextype jyend   = mysub*(jysub+1) - 1;

  /* Loop over the local array, computing the initial profile value.
     The global indices are (i,j) and the local indices are (iloc,jloc).
     Also set the id vector to zero for boundary points, one otherwise. */

  for (j = jybegin, jloc = 0; j <= jyend; j++, jloc++) {
    yfact = dy*j;
    for (i = ixbegin, iloc = 0; i <= ixend; i++, iloc++) {
      xfact = dx*i;
      loc = iloc + jloc*mxsub;
      u_1d(loc) = RCONST(16.0) * xfact * (ONE - xfact) * yfact * (ONE - yfact);

      if (i == 0 || i == data->mx - 1 || j == 0 || j == data->my - 1)
        id_1d(loc) = ZERO;
      else
        id_1d(loc) = ONE;
    }
  }

  /* Sync the device with the host */
  tu->sync<memory_space>();
  ti->sync<memory_space>();

  /* Initialize up. */

  N_VConst(ZERO, up);    /* Initially set up = 0. */

  /* resHeat sets res to negative of ODE RHS values at interior points. */
  resHeat(ZERO, uu, up, res, data);

  /* Copy -res into up to get correct initial up values on the device only! */
  N_VScale(-ONE, res, up);

  return(0);
}


/*
 * Print first lines of output and table heading
 */
static void PrintHeader(realtype rtol, realtype atol, UserData *data)
{
  printf("\nidaHeat2D_kry_p: Heat equation, parallel example problem for IDA\n");
  printf("            Discretized heat equation on 2D unit square.\n");
  printf("            Zero boundary conditions,");
  printf(" polynomial initial conditions.\n");
  printf("            Mesh dimensions: %d x %d", (int) data->mx, (int) data->my);
  printf("        Total system size: %ld\n\n", (long) data->mx * data->my);
  printf("Subgrid dimensions: %d x %d", (int) data->mxsub, (int) data->mysub);
  printf("        Processor array: %d x %d\n", (int) data->npex, (int) data->npey);
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("Tolerance parameters:  rtol = %Lg   atol = %Lg\n", rtol, atol);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("Tolerance parameters:  rtol = %g   atol = %g\n", rtol, atol);
#else
  printf("Tolerance parameters:  rtol = %g   atol = %g\n", rtol, atol);
#endif
  printf("Constraints set to force all solution components >= 0. \n");
  printf("SUPPRESSALG = SUNTRUE to suppress local error testing on ");
  printf("all boundary components. \n");
  printf("Linear solver: SUNSPGMR  ");
  printf("Preconditioner: diagonal elements only.\n");

  /* Print output table heading and initial line of table. */
  printf("\n   Output Summary (umax = max-norm of solution) \n\n");
  printf("  time     umax       k  nst  nni  nli   nre   nreLS    h      npe nps\n");
  printf("----------------------------------------------------------------------\n");
}


/*
 * PrintOutput: print max norm of solution and current solver statistics
 */
static void PrintOutput(int id, void *ida_mem, realtype t, N_Vector uu)
{
  realtype hused, umax;
  long int nst, nni, nje, nre, nreLS, nli, npe, nps;
  int kused, retval;

  umax = N_VMaxNorm(uu);

  if (id == 0) {

    retval = IDAGetLastOrder(ida_mem, &kused);
    check_retval(&retval, "IDAGetLastOrder", 1, id);
    retval = IDAGetNumSteps(ida_mem, &nst);
    check_retval(&retval, "IDAGetNumSteps", 1, id);
    retval = IDAGetNumNonlinSolvIters(ida_mem, &nni);
    check_retval(&retval, "IDAGetNumNonlinSolvIters", 1, id);
    retval = IDAGetNumResEvals(ida_mem, &nre);
    check_retval(&retval, "IDAGetNumResEvals", 1, id);
    retval = IDAGetLastStep(ida_mem, &hused);
    check_retval(&retval, "IDAGetLastStep", 1, id);
    retval = IDAGetNumJtimesEvals(ida_mem, &nje);
    check_retval(&retval, "IDAGetNumJtimesEvals", 1, id);
    retval = IDAGetNumLinIters(ida_mem, &nli);
    check_retval(&retval, "IDAGetNumLinIters", 1, id);
    retval = IDAGetNumLinResEvals(ida_mem, &nreLS);
    check_retval(&retval, "IDAGetNumLinResEvals", 1, id);
    retval = IDAGetNumPrecEvals(ida_mem, &npe);
    check_retval(&retval, "IDAGetPrecEvals", 1, id);
    retval = IDAGetNumPrecSolves(ida_mem, &nps);
    check_retval(&retval, "IDAGetNumPrecSolves", 1, id);

#if defined(SUNDIALS_EXTENDED_PRECISION)
    printf(" %5.2Lf %13.5Le  %d  %3ld  %3ld  %3ld  %4ld  %4ld  %9.2Le  %3ld %3ld\n",
           t, umax, kused, nst, nni, nje, nre, nreLS, hused, npe, nps);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    printf(" %5.2f %13.5e  %d  %3ld  %3ld  %3ld  %4ld  %4ld  %9.2e  %3ld %3ld\n",
           t, umax, kused, nst, nni, nje, nre, nreLS, hused, npe, nps);
#else
    printf(" %5.2f %13.5e  %d  %3ld  %3ld  %3ld  %4ld  %4ld  %9.2e  %3ld %3ld\n",
           t, umax, kused, nst, nni, nje, nre, nreLS, hused, npe, nps);
#endif

  }
}


/*
 * Print some final integrator statistics
 */
static void PrintFinalStats(void *ida_mem)
{
  long int netf, ncfn, ncfl;

  IDAGetNumErrTestFails(ida_mem, &netf);
  IDAGetNumNonlinSolvConvFails(ida_mem, &ncfn);
  IDAGetNumLinConvFails(ida_mem, &ncfl);

  printf("\nError test failures            = %ld\n", netf);
  printf("Nonlinear convergence failures = %ld\n", ncfn);
  printf("Linear convergence failures    = %ld\n", ncfl);
}

/*
 * Check function return value...
 *   opt == 0 means SUNDIALS function allocates memory so check if
 *            returned NULL pointer
 *   opt == 1 means SUNDIALS function returns a flag so check if
 *            flag >= 0
 *   opt == 2 means function allocates memory so check if returned
 *            NULL pointer
 */
static int check_retval(void *flagvalue, const char *funcname, int opt, int id)
{
  int *errflag;

  if (opt == 0 && flagvalue == NULL) {
    /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
    fprintf(stderr,
            "\nSUNDIALS_ERROR(%d): %s() failed - returned NULL pointer\n\n",
            id, funcname);
    return(1);
  } else if (opt == 1) {
    /* Check if flag < 0 */
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr,
              "\nSUNDIALS_ERROR(%d): %s() failed with flag = %d\n\n",
              id, funcname, *errflag);
      return(1);
    }
  } else if (opt == 2 && flagvalue == NULL) {
    /* Check if function returned NULL pointer - no memory allocated */
    fprintf(stderr,
            "\nMEMORY_ERROR(%d): %s() failed - returned NULL pointer\n\n",
            id, funcname);
    return(1);
  }

  return(0);
}
