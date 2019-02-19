/*
 * -----------------------------------------------------------------
 * Programmer(s): Slaven Peles @ LLNL
 * -----------------------------------------------------------------
 * Based on work by Daniel R. Reynolds @ SMU
 *         Allan Taylor, Alan Hindmarsh and Radu Serban @ LLNL
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
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <ida/ida.h>
#include <nvector/nvector_mpicuda.h>
#include <ida/ida_spils.h>
#include <sunlinsol/sunlinsol_spgmr.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_mpi_types.h>
#include <sundials/sundials_math.h>

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

typedef struct {
  int thispe, npex, npey, ixsub, jysub;
  sunindextype mx, my, mxsub, mysub;
  realtype     dx, dy, coeffx, coeffy, coeffxy;
  realtype    *uext; /* device array */
  realtype    *host_send_buff;
  realtype    *host_recv_buff;
  realtype    *dev_send_buff;
  realtype    *dev_recv_buff;
  N_Vector     pp;    /* vector of diagonal preconditioner elements */
  MPI_Comm  comm;
} *UserData;

/* User-supplied residual function and supporting routines */

int resHeat(realtype tt, N_Vector uu, N_Vector up,
            N_Vector rr, void *user_data);

static int rescomm(N_Vector uu, N_Vector up, void *user_data);

static int reslocal(realtype tt, N_Vector uu, N_Vector up,
                    N_Vector res,  void *user_data);

static int BSend(MPI_Comm comm, int thispe,
                 int ixsub, int jysub, int npex, int npey,
                 sunindextype mxsub, sunindextype mysub,
                 const realtype *uarray, realtype *dev_send_buff, realtype *host_send_buff);

static int BRecvPost(MPI_Comm comm, MPI_Request request[], int thispe,
                     int ixsub, int jysub, int npex, int npey,
                     sunindextype mxsub, sunindextype mysub,
                     realtype *host_recv_buff);

static int BRecvWait(MPI_Request request[],
                     int ixsub, int jysub, int npex, int npey,
                     sunindextype mxsub, sunindextype mysub,
                     realtype *uext, const realtype *host_recv_buff, realtype *dev_recv_buff);

/* User-supplied preconditioner routines */

int PsolveHeat(realtype tt, N_Vector uu, N_Vector up, N_Vector rr,
               N_Vector rvec, N_Vector zvec, realtype c_j,
               realtype delta, void *user_data);

int PsetupHeat(realtype tt, N_Vector yy, N_Vector yp, N_Vector rr,
               realtype c_j, void *user_data);

/* Private function to check function return values */

static int InitUserData(int thispe, MPI_Comm comm, UserData data);

static int AllocUserData(int thispe, MPI_Comm comm, N_Vector uu, UserData data);

static int DeleteUserData(UserData data);

static int SetInitialProfile(N_Vector uu, N_Vector up, N_Vector id,
                             N_Vector res, UserData data);

static void PrintHeader(realtype rtol, realtype atol, UserData data);

static void PrintOutput(int id, void *ida_mem, realtype t, N_Vector uu);

static void PrintFinalStats(void *ida_mem);

static int check_flag(void *flagvalue, const char *funcname, int opt, int id);


/*
 *--------------------------------------------------------------------
 * CUDA Kernels
 *--------------------------------------------------------------------
 */

__global__
void PsetupHeatKernel(realtype *ppv, sunindextype mx, sunindextype my,
                      sunindextype ibc, sunindextype jbc,
                      sunindextype i0, sunindextype j0, realtype pelinv)
{
  sunindextype i, j, tid;

  /* Loop over all grid points. */
  tid = blockDim.x * blockIdx.x + threadIdx.x;

  if (tid < (mx - ibc)*(my - jbc)) {
    i = tid % (mx - ibc) + i0;
    j = tid / (mx - ibc) + j0;

    /* Loop over interior points; ppv_i = 1/J_ii */
    ppv[i + j*mx] = pelinv;
  }
}


__global__
void CopyLocalToExtendedArray(const realtype *uuv, realtype *uext,
                              sunindextype mx, sunindextype my)
{
  sunindextype i, j, tid;

  /* Loop over all grid points. */
  tid = blockDim.x * blockIdx.x + threadIdx.x;

  if (tid < mx*my) {
    i = tid % mx;
    j = tid / mx;

    uext[(i+1) + (j+1)*(mx+2)] = uuv[i + j*mx];
  }
}


__global__
void LocalResidualKernel(const realtype *uext, const realtype *upv, realtype *resv,
                         sunindextype mx, sunindextype my,
                         sunindextype ibc, sunindextype jbc,
                         sunindextype i0, sunindextype j0,
                         realtype coeffx, realtype coeffy, realtype coeffxy)
{
  sunindextype i, j, tid;

  /* Loop over all grid points. */
  tid = blockDim.x * blockIdx.x + threadIdx.x;

  if (tid < (mx - ibc)*(my - jbc)) {
    i = tid % (mx - ibc) + i0;
    j = tid / (mx - ibc) + j0;

    sunindextype locu  = i + j*mx;
    sunindextype locue = (i+1) + (j+1)*(mx+2);

    realtype termx   = coeffx * (uext[locue-1]      + uext[locue+1]);
    realtype termy   = coeffy * (uext[locue-(mx+2)] + uext[locue+(mx+2)]);
    realtype termctr = coeffxy * uext[locue];
    resv[locu] = upv[locu] - (termx + termy - termctr);
  }
}


__global__
void CopyToBottomBuffer(const realtype *uarray, realtype *bufbottom,
                        sunindextype mx)
{
  sunindextype tid;

  /* Loop over all grid points. */
  tid = blockDim.x * blockIdx.x + threadIdx.x;

  if (tid < mx) {
      bufbottom[tid] = uarray[tid];
  }
}


__global__
void CopyToTopBuffer(const realtype *uarray, realtype *buftop,
                     sunindextype mx, sunindextype my)
{
  sunindextype tid;

  /* Loop over all grid points. */
  tid = blockDim.x * blockIdx.x + threadIdx.x;

  if (tid < mx) {
      buftop[tid] = uarray[(my-1)*mx + tid];
  }
}


__global__
void CopyToLeftBuffer(const realtype *uarray, realtype *bufleft,
                      sunindextype mx, sunindextype my)
{
  sunindextype tid;

  /* Loop over all grid points. */
  tid = blockDim.x * blockIdx.x + threadIdx.x;

  if (tid < my) {
      bufleft[tid] = uarray[tid*mx];
  }
}


__global__
void CopyToRightBuffer(const realtype *uarray, realtype *bufright,
                       sunindextype mx, sunindextype my)
{
  sunindextype tid;

  /* Loop over all grid points. */
  tid = blockDim.x * blockIdx.x + threadIdx.x;

  if (tid < my) {
      bufright[tid] = uarray[tid*mx + (mx-1)];
  }
}


__global__
void CopyFromBottomBuffer(const realtype *bufbottom, realtype *uext,
                          sunindextype mx)
{
  sunindextype tid;

  /* Loop over all grid points. */
  tid = blockDim.x * blockIdx.x + threadIdx.x;

  if (tid < mx) {
      uext[1 + tid] = bufbottom[tid];
  }
}


__global__
void CopyFromTopBuffer(const realtype *buftop, realtype *uext,
                       sunindextype mx, sunindextype my)
{
  sunindextype tid;

  /* Loop over all grid points. */
  tid = blockDim.x * blockIdx.x + threadIdx.x;

  if (tid < mx) {
      uext[(1 + (my+1)*(mx+2)) + tid] = buftop[tid];
  }
}


__global__
void CopyFromLeftBuffer(const realtype *bufleft, realtype *uext,
                        sunindextype mx, sunindextype my)
{
  sunindextype tid;

  /* Loop over all grid points. */
  tid = blockDim.x * blockIdx.x + threadIdx.x;

  if (tid < my) {
      uext[(tid+1)*(mx+2)] = bufleft[tid];
  }
}


__global__
void CopyFromRightBuffer(const realtype *bufright, realtype *uext,
                         sunindextype mx, sunindextype my)
{
  sunindextype tid;

  /* Loop over all grid points. */
  tid = blockDim.x * blockIdx.x + threadIdx.x;

  if (tid < my) {
      uext[(tid+2)*(mx+2) - 1] = bufright[tid];
  }
}


/*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 */

int main(int argc, char *argv[])
{
  MPI_Comm comm;
  void *ida_mem;
  SUNLinearSolver LS;
  UserData data;
  int iout, thispe, ier, npes;
  sunindextype Neq, local_N;
  realtype rtol, atol, t0, t1, tout, tret;
  N_Vector uu, up, constraints, id, res;

  ida_mem = NULL;
  LS = NULL;
  data = NULL;
  uu = up = constraints = id = res = NULL;

  /* Get processor number and total number of pe's. */

  MPI_Init(&argc, &argv);
  comm = MPI_COMM_WORLD;
  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &thispe);

  /* Allocate and initialize the data structure */
  data = (UserData) malloc(sizeof *data);
  if(check_flag((void *)data, "malloc", 2, thispe))
    MPI_Abort(comm, 1);

  InitUserData(thispe, comm, data);

  /* Check if the number of MPI processes matches the number of subgrids */
  if (npes != (data->npex * data->npey)) {
    if (thispe == 0)
      fprintf(stderr,
              "\nMPI_ERROR(0): npes = %d is not equal to NPEX*NPEY = %d\n",
              npes, data->npex * data->npey);
    free(data);
    MPI_Finalize();
    return(1);
  }

  /* Set local length local_N and global length Neq. */

  local_N = data->mxsub * data->mysub;
  Neq     = data->mx * data->my;

  /* Allocate and initialize N-vectors. */

  uu = N_VNew_Cuda(comm, local_N, Neq);
  if(check_flag((void *)uu, "N_VNew_Parallel", 0, thispe))
    MPI_Abort(comm, 1);

  up = N_VClone(uu);
  if(check_flag((void *)up, "N_VClone", 0, thispe))
    MPI_Abort(comm, 1);

  res = N_VClone(uu);
  if(check_flag((void *)res, "N_VClone", 0, thispe))
    MPI_Abort(comm, 1);

  constraints = N_VClone(uu);
  if(check_flag((void *)constraints, "N_VClone", 0, thispe))
    MPI_Abort(comm, 1);

  id = N_VClone(uu);
  if(check_flag((void *)id, "N_VClone", 0, thispe))
    MPI_Abort(comm, 1);

  /* Allocate user data extended vector and MPI buffers */
  ier = AllocUserData(thispe, comm, uu, data);
  if(check_flag(&ier, "AllocUserData", 1, thispe)) MPI_Abort(comm, 1);


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
  if(check_flag((void *)ida_mem, "IDACreate", 0, thispe)) MPI_Abort(comm, 1);

  ier = IDASetUserData(ida_mem, data);
  if(check_flag(&ier, "IDASetUserData", 1, thispe)) MPI_Abort(comm, 1);

  ier = IDASetSuppressAlg(ida_mem, SUNTRUE);
  if(check_flag(&ier, "IDASetSuppressAlg", 1, thispe)) MPI_Abort(comm, 1);

  ier = IDASetId(ida_mem, id);
  if(check_flag(&ier, "IDASetId", 1, thispe)) MPI_Abort(comm, 1);

  ier = IDASetConstraints(ida_mem, constraints);
  if(check_flag(&ier, "IDASetConstraints", 1, thispe)) MPI_Abort(comm, 1);
  N_VDestroy(constraints);

  ier = IDAInit(ida_mem, resHeat, t0, uu, up);
  if(check_flag(&ier, "IDAInit", 1, thispe)) MPI_Abort(comm, 1);

  ier = IDASStolerances(ida_mem, rtol, atol);
  if(check_flag(&ier, "IDASStolerances", 1, thispe)) MPI_Abort(comm, 1);

  /* Call SUNSPGMR and IDASetLinearSolver to specify the linear solver. */

  LS = SUNSPGMR(uu, PREC_LEFT, 0);  /* use default maxl */
  if(check_flag((void *)LS, "SUNSPGMR", 0, thispe)) MPI_Abort(comm, 1);

  ier = IDASpilsSetLinearSolver(ida_mem, LS);
  if(check_flag(&ier, "IDASpilsSetLinearSolver", 1, thispe)) MPI_Abort(comm, 1);

  ier = IDASpilsSetPreconditioner(ida_mem, PsetupHeat, PsolveHeat);
  if(check_flag(&ier, "IDASpilsSetPreconditioner", 1, thispe)) MPI_Abort(comm, 1);

  /* Print output heading (on processor 0 only) and intial solution  */

  if (thispe == 0) PrintHeader(rtol, atol, data);
  PrintOutput(thispe, ida_mem, t0, uu);

  /* Loop over tout, call IDASolve, print output. */

  for (tout = t1, iout = 1; iout <= NOUT; iout++, tout *= TWO) {

    ier = IDASolve(ida_mem, tout, &tret, uu, up, IDA_NORMAL);
    if(check_flag(&ier, "IDASolve", 1, thispe)) MPI_Abort(comm, 1);

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

  DeleteUserData(data);
  free(data);

  MPI_Finalize();

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
  sunindextype ibc, i0, jbc, j0;

  /* Unwrap the user data */
  UserData data = (UserData) user_data;
  const int ixsub = data->ixsub;
  const int jysub = data->jysub;
  const int npex  = data->npex;
  const int npey  = data->npey;
  const sunindextype mxsub = data->mxsub;
  const sunindextype mysub = data->mysub;
  realtype *ppv = N_VGetDeviceArrayPointer_Cuda(data->pp);

  /* Calculate the value for the inverse element of the diagonal preconditioner */
  const realtype pelinv = ONE/(c_j + data->coeffxy);

  /* Initially set all pp elements on the device to one. */
  N_VConst(ONE, data->pp);

  ibc = (ixsub == 0) || (ixsub == npex-1);
  i0  = (ixsub == 0);
  jbc = (jysub == 0) || (jysub == npey-1);
  j0  = (jysub == 0);

  unsigned block = 256;
  unsigned grid = ((mxsub - ibc)*(mysub - jbc) + block - 1) / block;

  PsetupHeatKernel<<<grid, block>>>(ppv, mxsub, mysub, ibc, jbc, i0, j0, pelinv);

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
  UserData data = (UserData) user_data;

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
  UserData data = (UserData) user_data;

  /* Get comm, thispe, subgrid indices, data sizes */
  MPI_Comm comm = data->comm;
  const int thispe = data->thispe;
  const int ixsub  = data->ixsub;
  const int jysub  = data->jysub;
  const int npex   = data->npex;
  const int npey   = data->npey;
  const sunindextype mxsub = data->mxsub;
  const sunindextype mysub = data->mysub;

  /* Get pointers to buffers and extended solution vector data array uext. */
  realtype *uext = data->uext;
  realtype *host_send_buff = data->host_send_buff;
  realtype *host_recv_buff = data->host_recv_buff;
  realtype *dev_send_buff  = data->dev_send_buff;
  realtype *dev_recv_buff  = data->dev_recv_buff;

  /* Get solution vector data. */
  const realtype *uarray = N_VGetDeviceArrayPointer_Cuda(uu);

  /* Set array of MPI requests */
  MPI_Request request[4];

  /* Start receiving boundary data from neighboring PEs. */
  BRecvPost(comm, request, thispe, ixsub, jysub, npex, npey, mxsub, mysub, host_recv_buff);

  /* Send data from boundary of local grid to neighboring PEs. */
  BSend(comm, thispe, ixsub, jysub, npex, npey, mxsub, mysub, uarray, dev_send_buff, host_send_buff);

  /* Finish receiving boundary data from neighboring PEs. */
  BRecvWait(request, ixsub, jysub, npex, npey, mxsub, mysub, uext, host_recv_buff, dev_recv_buff);

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
  UserData data = (UserData) user_data;

  /* Get subgrid indices, array sizes, and grid coefficients. */
  const int ixsub = data->ixsub;
  const int jysub = data->jysub;
  const int npex  = data->npex;
  const int npey  = data->npey;
  const sunindextype mxsub  = data->mxsub;
  const sunindextype mysub  = data->mysub;
  const realtype coeffx  = data->coeffx;
  const realtype coeffy  = data->coeffy;
  const realtype coeffxy = data->coeffxy;

  /* Vector data arrays, extended work array uext. */
  const realtype *uuv = N_VGetDeviceArrayPointer_Cuda(uu);
  const realtype *upv = N_VGetDeviceArrayPointer_Cuda(up);
  realtype *resv = N_VGetDeviceArrayPointer_Cuda(rr);
  realtype *uext = data->uext;

  sunindextype ibc, i0, jbc, j0;

  /* Initialize all elements of rr to uu. This sets the boundary
     elements simply without indexing hassles. */

  N_VScale(ONE, uu, rr);

  /* Copy local segment of u vector into the working extended array uext.
     This completes uext prior to the computation of the rr vector.
     uext and uuv must be on the device.     */
  unsigned block = 256;
  unsigned grid = (mxsub*mysub + block - 1) / block;

  CopyLocalToExtendedArray<<<grid, block>>>(uuv, uext, mxsub, mysub);

  /* Set loop limits for the interior of the local subgrid. */
  ibc = (ixsub == 0) || (ixsub == npex-1);
  i0  = (ixsub == 0);
  jbc = (jysub == 0) || (jysub == npey-1);
  j0  = (jysub == 0);

  /* Compute local residual; uext, upv, and resv must be on the device */
  block = 256;
  grid = ((mxsub - ibc)*(mysub - jbc) + block - 1) / block;

  LocalResidualKernel<<<grid, block>>>(uext, upv, resv, mxsub, mysub, ibc, jbc,
                                       i0, j0, coeffx, coeffy, coeffxy);

  return(0);

}

/*
 * Routine to send boundary data to neighboring PEs.
 */

static int BSend(MPI_Comm comm, int thispe,
                 int ixsub, int jysub, int npex, int npey,
                 sunindextype mxsub, sunindextype mysub,
                 const realtype *uarray, realtype *dev_send_buff, realtype *host_send_buff)
{
  cudaError_t err;
  //const sunindextype zero = 0;
  /* Have left, right, top and bottom device buffers use the same dev_send_buff. */
  realtype *d_bufleft   = dev_send_buff;
  realtype *d_bufright  = dev_send_buff + mysub;
  realtype *d_buftop    = dev_send_buff + 2*mysub;
  realtype *d_bufbottom = dev_send_buff + 2*mysub + mxsub;

  /* Have left, right, top and bottom host buffers use the same host_send_buff. */
  realtype *h_bufleft   = host_send_buff;
  realtype *h_bufright  = host_send_buff + mysub;
  realtype *h_buftop    = host_send_buff + 2*mysub;
  realtype *h_bufbottom = host_send_buff + 2*mysub + mxsub;

  /* If jysub > 0, send data from bottom x-line of u.  (via bufbottom) */

  if (jysub != 0) {
    // Device kernel here to copy from uarray to the buffer on the device
    unsigned block = 256;
    unsigned grid = (mxsub + block - 1) / block;
    CopyToBottomBuffer<<<grid, block>>>(uarray, d_bufbottom, mxsub);

    // Copy buffer to the host
    err = cudaMemcpy(h_bufbottom, d_bufbottom, mxsub*sizeof(realtype), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
      printf("Bottom buffer: Copy from device to host failed with code %d... \n", err);
      printf("%ld %ld\n", h_bufbottom, d_bufbottom);
      return -1;
    }
    // MPI send buffer
    MPI_Send(h_bufbottom, mxsub, PVEC_REAL_MPI_TYPE, thispe-npex, 0, comm);
  }

  /* If jysub < NPEY-1, send data from top x-line of u. (via buftop) */

  if (jysub != npey-1) {
    // Device kernel here to copy from uarray to the buffer on the device
    unsigned block = 256;
    unsigned grid = (mxsub + block - 1) / block;
    CopyToTopBuffer<<<grid, block>>>(uarray, d_buftop, mxsub, mysub);

    // Copy buffer to the host
    err = cudaMemcpy(h_buftop, d_buftop, mxsub*sizeof(realtype), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
      printf("Top buffer: Copy from device to host failed ... \n");
      return -1;
    }
    // MPI send buffer
    MPI_Send(h_buftop, mxsub, PVEC_REAL_MPI_TYPE, thispe+npex, 0, comm);
  }

  /* If ixsub > 0, send data from left y-line of u (via bufleft). */

  if (ixsub != 0) {
    // Device kernel here to copy from uarray to the buffer on the device
    unsigned block = 256;
    unsigned grid = (mysub + block - 1) / block;
    CopyToLeftBuffer<<<grid, block>>>(uarray, d_bufleft, mxsub, mysub);

    // Copy buffer to the host
    err = cudaMemcpy(h_bufleft, d_bufleft, mysub*sizeof(realtype), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
      printf("Left buffer: Copy from device to host failed ... \n");
      return -1;
    }
    // MPI send buffer
    MPI_Send(h_bufleft, mysub, PVEC_REAL_MPI_TYPE, thispe-1, 0, comm);
  }

  /* If ixsub < NPEX-1, send data from right y-line of u (via bufright). */

  if (ixsub != npex-1) {
    // Device kernel here to copy from uarray to the buffer on the device
    unsigned block = 256;
    unsigned grid = (mysub + block - 1) / block;
    CopyToRightBuffer<<<grid, block>>>(uarray, d_bufright, mxsub, mysub);

    // Copy buffer to the host
    err = cudaMemcpy(h_bufright, d_bufright, mysub*sizeof(realtype), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
      printf("Right buffer: Copy from device to host failed ... \n");
      return -1;
    }
    // MPI send buffer
    MPI_Send(h_bufright, mysub, PVEC_REAL_MPI_TYPE, thispe+1, 0, comm);
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

static int BRecvPost(MPI_Comm comm, MPI_Request request[], int thispe,
                     int ixsub, int jysub, int npex, int npey,
                     sunindextype mxsub, sunindextype mysub,
                     realtype *host_recv_buff)
{
  /* Have left, right, top and bottom buffers use the same host_recv_buff. */
  realtype *bufleft   = host_recv_buff;
  realtype *bufright  = host_recv_buff + mysub;
  realtype *buftop    = host_recv_buff + 2*mysub;
  realtype *bufbottom = host_recv_buff + 2*mysub + mxsub;

  /* If jysub > 0, receive data for bottom x-line of uext. */
  if (jysub != 0) {
    MPI_Irecv(bufbottom, mxsub, PVEC_REAL_MPI_TYPE,
              thispe-npex, 0, comm, &request[0]);
  }

  /* If jysub < NPEY-1, receive data for top x-line of uext. */
  if (jysub != npey-1) {
    MPI_Irecv(buftop, mxsub, PVEC_REAL_MPI_TYPE,
              thispe+npex, 0, comm, &request[1]);
  }

  /* If ixsub > 0, receive data for left y-line of uext (via bufleft). */
  if (ixsub != 0) {
    MPI_Irecv(&bufleft[0], mysub, PVEC_REAL_MPI_TYPE,
              thispe-1, 0, comm, &request[2]);
  }

  /* If ixsub < NPEX-1, receive data for right y-line of uext (via bufright). */
  if (ixsub != npex-1) {
    MPI_Irecv(&bufright[0], mysub, PVEC_REAL_MPI_TYPE,
              thispe+1, 0, comm, &request[3]);
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

static int BRecvWait(MPI_Request request[], int ixsub, int jysub,
                     int npex, int npey,
                     sunindextype mxsub, sunindextype mysub,
                     realtype *uext, const realtype *host_recv_buff, realtype *dev_recv_buff)
{
  cudaError_t err;
  MPI_Status status;

  const realtype *h_bufleft   = host_recv_buff;
  const realtype *h_bufright  = host_recv_buff + mysub;
  const realtype *h_buftop    = host_recv_buff + 2*mysub;
  const realtype *h_bufbottom = host_recv_buff + 2*mysub + mxsub;

  realtype *d_bufleft   = dev_recv_buff;
  realtype *d_bufright  = dev_recv_buff + mysub;
  realtype *d_buftop    = dev_recv_buff + 2*mysub;
  realtype *d_bufbottom = dev_recv_buff + 2*mysub + mxsub;

  /* If jysub > 0, receive data for bottom x-line of uext. */
  if (jysub != 0) {
    MPI_Wait(&request[0], &status);
    /* Copy the buffer from the host to the device */
    err = cudaMemcpy(d_bufbottom, h_bufbottom, mxsub*sizeof(realtype), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
      printf("Copy from host to device failed ... \n");
      return -1;
    }
    /* Copy the bottom dev_recv_buff to uext. */
    unsigned block = 256;
    unsigned grid = (mxsub + block - 1) / block;
    CopyFromBottomBuffer<<<grid, block>>>(d_bufbottom, uext, mxsub);
  }

  /* If jysub < NPEY-1, receive data for top x-line of uext. */
  if (jysub != npey-1) {
    MPI_Wait(&request[1], &status);
    /* Copy the buffer from the host to the device */
    err = cudaMemcpy(d_buftop, h_buftop, mxsub*sizeof(realtype), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
      printf("Copy from host to device failed ... \n");
      return -1;
    }
    /* Copy the top dev_recv_buff to uext. */
    unsigned block = 256;
    unsigned grid = (mxsub + block - 1) / block;
    CopyFromTopBuffer<<<grid, block>>>(d_buftop, uext, mxsub, mysub);
  }

  /* If ixsub > 0, receive data for left y-line of uext (via bufleft). */
  if (ixsub != 0) {
    MPI_Wait(&request[2], &status);
    /* Copy the buffer from the host to the device */
    err = cudaMemcpy(d_bufleft, h_bufleft, mysub*sizeof(realtype), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
      printf("Copy from host to device failed ... \n");
      return -1;
    }
    /* Copy the left dev_recv_buff to uext. */
    unsigned block = 256;
    unsigned grid = (mysub + block - 1) / block;
    CopyFromLeftBuffer<<<grid, block>>>(d_bufleft, uext, mxsub, mysub);
  }

  /* If ixsub < NPEX-1, receive data for right y-line of uext (via bufright). */
  if (ixsub != npex-1) {
    MPI_Wait(&request[3], &status);
    /* Copy the buffer from the host to the device */
    err = cudaMemcpy(d_bufright, h_bufright, mysub*sizeof(realtype), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
      printf("Copy from host to device failed ... \n");
      return -1;
    }
    /* Copy the right dev_recv_buff to uext. */
    unsigned block = 256;
    unsigned grid = (mysub + block - 1) / block;
    CopyFromRightBuffer<<<grid, block>>>(d_bufright, uext, mxsub, mysub);
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

static int InitUserData(int thispe, MPI_Comm comm, UserData data)
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

  data->uext = NULL;
  data->host_send_buff = NULL;
  data->host_recv_buff = NULL;
  data->dev_send_buff  = NULL;
  data->dev_recv_buff  = NULL;

  return(0);
}


/*
 * AllocUserData allocates memory for the extended vector uext
 * and MPI communication buffers.
 */

static int AllocUserData(int thispe, MPI_Comm comm, N_Vector uu, UserData data)
{
  cudaError_t err;
  sunindextype mxsub = data->mxsub;
  sunindextype mysub = data->mysub;

  /* An N-vector to hold preconditioner. */
  data->pp = N_VClone(uu);
  if(data->pp == NULL) {
    MPI_Abort(comm, 1);
    return -1;
  }

  /* Allocate local extended vector (includes ghost nodes) */
  err = cudaMalloc((void**) &data->uext, (mxsub + 2)*(mysub +2)*sizeof(realtype));
  if(err != cudaSuccess) {
    printf("Failed to allocate uext ... \n");
    N_VDestroy(data->pp);
    MPI_Abort(comm, 1);
    return -1;
  }

  /* Allocate local host send buffer */
  data->host_send_buff = (realtype*) malloc(2*(mxsub + mysub)*sizeof(realtype));
  if(data->host_send_buff == NULL) {
    N_VDestroy(data->pp);
    free(data->uext);
    MPI_Abort(comm, 1);
    return -1;
  }

  data->host_recv_buff = (realtype*) malloc(2*(mxsub + mysub)*sizeof(realtype));
  if(data->host_recv_buff == NULL) {
    N_VDestroy(data->pp);
    free(data->uext);
    free(data->host_send_buff);
    MPI_Abort(comm, 1);
    return -1;
  }

  /* Allocate local device send buffer */
  err = cudaMalloc((void**) &data->dev_send_buff, 2*(mxsub + mysub)*sizeof(realtype));
  if(err != cudaSuccess) {
    printf("Failed to allocate dev_send_buff ... \n");
    N_VDestroy(data->pp);
    cudaFree(data->uext);
    free(data->host_send_buff);
    free(data->host_recv_buff);
    MPI_Abort(comm, 1);
    return -1;
  }

  /* Allocate local device send buffer */
  err = cudaMalloc((void**) &data->dev_recv_buff, 2*(mxsub + mysub)*sizeof(realtype));
  if(err != cudaSuccess) {
    printf("Failed to allocate dev_recv_buff ... \n");
    N_VDestroy(data->pp);
    cudaFree(data->uext);
    free(data->host_send_buff);
    free(data->host_recv_buff);
    cudaFree(data->dev_send_buff);
    MPI_Abort(comm, 1);
    return -1;
  }

  return 0;
}


static int DeleteUserData(UserData data)
{
  if (data->pp != NULL)
    N_VDestroy(data->pp);
  if (data->uext != NULL)
    cudaFree(data->uext);
  if (data->host_send_buff != NULL)
    free(data->host_send_buff);
  if (data->host_recv_buff != NULL)
    free(data->host_recv_buff);
  if (data->dev_send_buff != NULL)
    cudaFree(data->dev_send_buff);
  if (data->dev_recv_buff != NULL)
    cudaFree(data->dev_recv_buff);
  return 0;
}

/*
 * SetInitialProfile sets the initial values for the problem.
 */

static int SetInitialProfile(N_Vector uu, N_Vector up,  N_Vector id,
                             N_Vector res, UserData data)
{
  sunindextype i, iloc, j, jloc, loc;
  realtype xfact, yfact;

  /* Initialize uu. */

  // Get host pointer
  realtype *uudata = N_VGetHostArrayPointer_Cuda(uu);
  realtype *iddata = N_VGetHostArrayPointer_Cuda(id);

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
      uudata[loc] = RCONST(16.0) * xfact * (ONE - xfact) * yfact * (ONE - yfact);

      if (i == 0 || i == data->mx - 1 || j == 0 || j == data->my - 1)
        iddata[loc] = ZERO;
      else
        iddata[loc] = ONE;
    }
  }

  // Synchronize data from the host to the device for uu and id vectors
  N_VCopyToDevice_Cuda(uu);
  N_VCopyToDevice_Cuda(id);

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

static void PrintHeader(realtype rtol, realtype atol, UserData data)
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
  int kused, ier;

  umax = N_VMaxNorm(uu);

  if (id == 0) {

    ier = IDAGetLastOrder(ida_mem, &kused);
    check_flag(&ier, "IDAGetLastOrder", 1, id);
    ier = IDAGetNumSteps(ida_mem, &nst);
    check_flag(&ier, "IDAGetNumSteps", 1, id);
    ier = IDAGetNumNonlinSolvIters(ida_mem, &nni);
    check_flag(&ier, "IDAGetNumNonlinSolvIters", 1, id);
    ier = IDAGetNumResEvals(ida_mem, &nre);
    check_flag(&ier, "IDAGetNumResEvals", 1, id);
    ier = IDAGetLastStep(ida_mem, &hused);
    check_flag(&ier, "IDAGetLastStep", 1, id);
    ier = IDASpilsGetNumJtimesEvals(ida_mem, &nje);
    check_flag(&ier, "IDASpilsGetNumJtimesEvals", 1, id);
    ier = IDASpilsGetNumLinIters(ida_mem, &nli);
    check_flag(&ier, "IDASpilsGetNumLinIters", 1, id);
    ier = IDASpilsGetNumResEvals(ida_mem, &nreLS);
    check_flag(&ier, "IDASpilsGetNumResEvals", 1, id);
    ier = IDASpilsGetNumPrecEvals(ida_mem, &npe);
    check_flag(&ier, "IDASpilsGetPrecEvals", 1, id);
    ier = IDASpilsGetNumPrecSolves(ida_mem, &nps);
    check_flag(&ier, "IDASpilsGetNumPrecSolves", 1, id);

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
  IDASpilsGetNumConvFails(ida_mem, &ncfl);

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

static int check_flag(void *flagvalue, const char *funcname, int opt, int id)
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
