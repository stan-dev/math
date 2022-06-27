/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * Shared implementaiton file for 2D diffusion benchmark problem
 * ---------------------------------------------------------------------------*/

#include "diffusion_2D.hpp"

// -----------------------------------------------------------------------------
// ODE and DAE problem defining functions
// -----------------------------------------------------------------------------

#if defined(BENCHMARK_ODE)

int diffusion(realtype t, N_Vector u, N_Vector f, void *user_data)
{
#ifdef SUNDIALS_BUILD_WITH_PROFILING
  // Access problem data
  UserData *udata = (UserData *) user_data;
#endif

  SUNDIALS_CXX_MARK_FUNCTION(udata->prof);

  // Compute the Laplacian
  int flag = laplacian(t, u, f, user_data);
  if (check_flag(&flag, "laplacian", 1))
    return -1;


  return 0;
}

#elif defined(BENCHMARK_DAE)

int diffusion(realtype t, N_Vector u, N_Vector up, N_Vector res,
              void *user_data)
{
#ifdef SUNDIALS_BUILD_WITH_PROFILING
  // Access problem data
  UserData *udata = (UserData *) user_data;
#endif

  SUNDIALS_CXX_MARK_FUNCTION(udata->prof);

  // Compute the Laplacian
  int flag = laplacian(t, u, res, user_data);
  if (check_flag(&flag, "laplacian", 1))
    return -1;

  // Compute the residual
  N_VLinearSum(ONE, up, -ONE, res, res);


  return 0;
}

#else
#error "Missing ODE/DAE preprocessor directive"
#endif

// -----------------------------------------------------------------------------
// UserData public functions
// -----------------------------------------------------------------------------

// Parse command line inputs
int UserData::parse_args(vector<string> &args, bool outproc)
{
  vector<string>::iterator it;

  it = find(args.begin(), args.end(), "--help");
  if (it != args.end())
  {
    if (outproc) help();
    return 0;
  }

  it = find(args.begin(), args.end(), "--npx");
  if (it != args.end())
  {
    npx = stoi(*(it + 1));
    args.erase(it, it + 2);
  }

  it = find(args.begin(), args.end(), "--npy");
  if (it != args.end())
  {
    npy = stoi(*(it + 1));
    args.erase(it, it + 2);
  }

  it = find(args.begin(), args.end(), "--nx");
  if (it != args.end())
  {
    nx = stoi(*(it + 1));
    args.erase(it, it + 2);
  }

  it = find(args.begin(), args.end(), "--ny");
  if (it != args.end())
  {
    ny = stoi(*(it + 1));
    args.erase(it, it + 2);
  }

  it = find(args.begin(), args.end(), "--xu");
  if (it != args.end())
  {
    xu = stod(*(it + 1));
    args.erase(it, it + 2);
  }

  it = find(args.begin(), args.end(), "--yu");
  if (it != args.end())
  {
    yu = stod(*(it + 1));
    args.erase(it, it + 2);
  }

  it = find(args.begin(), args.end(), "--kx");
  if (it != args.end())
  {
    kx = stod(*(it + 1));
    args.erase(it, it + 2);
  }

  it = find(args.begin(), args.end(), "--ky");
  if (it != args.end())
  {
    ky = stod(*(it + 1));
    args.erase(it, it + 2);
  }

  it = find(args.begin(), args.end(), "--tf");
  if (it != args.end())
  {
    tf = stod(*(it + 1));
    args.erase(it, it + 2);
  }

  it = find(args.begin(), args.end(), "--noforcing");
  if (it != args.end())
  {
    forcing = false;
    args.erase(it);
  }

  // Recompute total number of nodes
  nodes = nx * ny;

  // Recompute x and y mesh spacing
  dx = xu / (nx - 1);
  dy = yu / (ny - 1);

  return 0;
}


// Print command line options
void UserData::help()
{
  cout << endl;
  cout << "Problem setup command line options:" << endl;
  cout << "  --nx <nx>    : x-direction mesh points" << endl;
  cout << "  --ny <ny>    : y-direction mesh points" << endl;
  cout << "  --xu <xu>    : x-direction upper bound" << endl;
  cout << "  --yu <yu>    : y-direction upper bound" << endl;
  cout << "  --kx <kx>    : x-direction diffusion coefficient" << endl;
  cout << "  --ky <kx>    : y-direction diffusion coefficient" << endl;
  cout << "  --noforcing  : disable forcing term" << endl;
  cout << "  --tf <time>  : final time" << endl;
}


void UserData::print()
{
  cout << endl;
  cout << " Problem options:" << endl;
  cout << " --------------------------------- " << endl;
  cout << "  nprocs         = " << np       << endl;
  cout << "  npx            = " << npx      << endl;
  cout << "  npy            = " << npy      << endl;
  cout << " --------------------------------- " << endl;
  cout << "  kx             = " << kx      << endl;
  cout << "  ky             = " << ky      << endl;
  cout << "  forcing        = " << forcing << endl;
  cout << "  tf             = " << tf      << endl;
  cout << "  xu             = " << xu      << endl;
  cout << "  yu             = " << yu      << endl;
  cout << "  nx             = " << nx      << endl;
  cout << "  ny             = " << ny      << endl;
  cout << "  nxl (proc 0)   = " << nx_loc  << endl;
  cout << "  nyl (proc 0)   = " << ny_loc  << endl;
  cout << "  dx             = " << dx      << endl;
  cout << "  dy             = " << dy      << endl;
  cout << " --------------------------------- " << endl;
}

int UserData::setup()
{
  int flag;

  // Check that this has not been called before
  if (Erecv != NULL || Wrecv != NULL ||
      Srecv != NULL || Nrecv != NULL)
  {
    cerr << "SetupDecomp error: parallel decomposition already set up" << endl;
    return -1;
  }

  // Get the number of processes
  flag = MPI_Comm_size(MPI_COMM_WORLD, &np);
  if (flag != MPI_SUCCESS)
  {
    cerr << "Error in MPI_Comm_size = " << flag << endl;
    return -1;
  }

  // Set up 2D Cartesian communicator
  int dims[2];
  dims[0] = (npx > 0) ? npx : 0;
  dims[1] = (npy > 0) ? npy : 0;

  int periods[2];
  periods[0] = 0;
  periods[1] = 0;

  flag = MPI_Dims_create(np, 2, dims);
  if (flag != MPI_SUCCESS)
  {
    cerr << "Error in MPI_Dims_create = " << flag << endl;
    return -1;
  }

  npx = dims[0];
  npy = dims[1];

  flag = MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &comm_c);
  if (flag != MPI_SUCCESS)
  {
    cerr << "Error in MPI_Cart_create = " << flag << endl;
    return -1;
  }

  // Get my rank in the new Cartesian communicator
  flag = MPI_Comm_rank(comm_c, &myid_c);
  if (flag != MPI_SUCCESS)
  {
    cerr << "Error in MPI_Comm_rank = " << flag << endl;
    return -1;
  }

  // Get dimension of the Cartesian communicator and my coordinates
  int coords[2];
  flag = MPI_Cart_get(comm_c, 2, dims, periods, coords);
  if (flag != MPI_SUCCESS)
  {
    cerr << "Error in MPI_Cart_get = " << flag << endl;
    return -1;
  }

  // Determine local extents in x-direction
  int idx         = coords[0];
  sunindextype qx = nx / dims[0];
  sunindextype rx = nx % dims[0];

  is = qx * idx + (idx < rx ? idx : rx);
  ie = is + qx - 1 + (idx < rx ? 1 : 0);

  // Sanity check
  if (ie > (nx - 1))
  {
    cerr << "Error ie > nx - 1" << endl;
    return -1;
  }

  // Determine local extents in y-direction
  int idy         = coords[1];
  sunindextype qy = ny / dims[1];
  sunindextype ry = ny % dims[1];

  js = qy * idy + (idy < ry ? idy : ry);
  je = js + qy - 1 + (idy < ry ? 1 : 0);

  // Sanity check
  if (je > (ny - 1))
  {
    cerr << "Error je > ny - 1" << endl;
    return -1;
  }

  // Number of local nodes
  nx_loc = (ie) - (is) + 1;
  ny_loc = (je) - (js) + 1;

  // Initialize global and local vector lengths
  nodes     = nx * ny;
  nodes_loc = nx_loc * ny_loc;

  // Determine if this proc has neighbors
  HaveNbrW = (is != 0);
  HaveNbrE = (ie != nx-1);
  HaveNbrS = (js != 0);
  HaveNbrN = (je != ny-1);

  // Allocate exchange buffers if necessary
  flag = allocate_buffers();
  if (flag)
  {
    cerr << "Error in AlocateBuffers = " << flag << endl;
    return -1;
  }

  // MPI neighborhood information
  int nbcoords[2];

  // West neighbor
  if (HaveNbrW)
  {
    nbcoords[0] = coords[0]-1;
    nbcoords[1] = coords[1];
    flag = MPI_Cart_rank(comm_c, nbcoords, &ipW);
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Cart_rank = " << flag << endl;
      return -1;
    }
  }

  // East neighbor
  if (HaveNbrE)
  {
    nbcoords[0] = coords[0]+1;
    nbcoords[1] = coords[1];
    flag = MPI_Cart_rank(comm_c, nbcoords, &ipE);
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Cart_rank = " << flag << endl;
      return -1;
    }
  }

  // South neighbor
  if (HaveNbrS)
  {
    nbcoords[0] = coords[0];
    nbcoords[1] = coords[1]-1;
    flag = MPI_Cart_rank(comm_c, nbcoords, &ipS);
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Cart_rank = " << flag << endl;
      return -1;
    }
  }

  // North neighbor
  if (HaveNbrN)
  {
    nbcoords[0] = coords[0];
    nbcoords[1] = coords[1]+1;
    flag = MPI_Cart_rank(comm_c, nbcoords, &ipN);
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Cart_rank = " << flag << endl;
      return -1;
    }
  }

  // Return success
  return 0;
}


// -----------------------------------------------------------------------------
// UserData boundary exchange functions
// -----------------------------------------------------------------------------


int UserData::start_exchange(const N_Vector u)
{
  int flag;

  SUNDIALS_CXX_MARK_FUNCTION(prof);

  // -------------
  // Pack buffers
  // -------------

  flag = pack_buffers(u);
  if (flag)
  {
    cerr << "Error in PackBuffers = " << flag << endl;
    return -1;
  }

  // -----------
  // Post Irecv
  // -----------

  if (HaveNbrW)
  {
    flag = MPI_Irecv(Wrecv, (int) ny_loc, MPI_SUNREALTYPE,
                     ipW, MPI_ANY_TAG, comm_c, &reqRW);
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Irecv = " << flag << endl;
      return -1;
    }
  }

  if (HaveNbrE)
  {
    flag = MPI_Irecv(Erecv, (int) ny_loc, MPI_SUNREALTYPE,
                     ipE, MPI_ANY_TAG, comm_c, &reqRE);
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Irecv = " << flag << endl;
      return -1;
    }
  }

  if (HaveNbrS)
  {
    flag = MPI_Irecv(Srecv, (int) nx_loc, MPI_SUNREALTYPE,
                     ipS, MPI_ANY_TAG, comm_c, &reqRS);
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Irecv = " << flag << endl;
      return -1;
    }
  }

  if (HaveNbrN)
  {
    flag = MPI_Irecv(Nrecv, (int) nx_loc, MPI_SUNREALTYPE,
                     ipN, MPI_ANY_TAG, comm_c, &reqRN);
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Irecv = " << flag << endl;
      return -1;
    }
  }

  // ----------
  // Send data
  // ----------

  // ensure packing has finished
  flag = DeviceSynchronize();
  if (flag != 0)
  {
    cerr << "Error in DeviceSynchronize" << endl;
    return -1;
  }

  if (HaveNbrW)
  {
    flag = MPI_Isend(Wsend, (int) ny_loc, MPI_SUNREALTYPE,
                     ipW, 0, comm_c, &reqSW);
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Isend = " << flag << endl;
      return -1;
    }
  }

  if (HaveNbrE)
  {
    flag = MPI_Isend(Esend, (int) ny_loc, MPI_SUNREALTYPE,
                     ipE, 1, comm_c, &reqSE);
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Isend = " << flag << endl;
      return -1;
    }
  }

  if (HaveNbrS)
  {
    flag = MPI_Isend(Ssend, (int) nx_loc, MPI_SUNREALTYPE,
                     ipS, 2, comm_c, &reqSS);
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Isend = " << flag << endl;
      return -1;
    }
  }

  if (HaveNbrN)
  {
    flag = MPI_Isend(Nsend, (int) nx_loc, MPI_SUNREALTYPE,
                     ipN, 3, comm_c, &reqSN);
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Isend = " << flag << endl;
      return -1;
    }
  }

  // Return success
  return 0;
}


int UserData::end_exchange()
{
  // Local variables
  int flag;
  MPI_Status stat;

  SUNDIALS_CXX_MARK_FUNCTION(prof);

  // Wait for messages to finish
  if (HaveNbrW)
  {
    flag = MPI_Wait(&reqRW, &stat);
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Wait = " << flag << endl;
      return -1;
    }
    flag = MPI_Wait(&reqSW, &stat);
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Wait = " << flag << endl;
      return -1;
    }
  }

  if (HaveNbrE)
  {
    flag = MPI_Wait(&reqRE, &stat);
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Wait = " << flag << endl;
      return -1;
    }
    flag = MPI_Wait(&reqSE, &stat);
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Wait = " << flag << endl;
      return -1;
    }
  }

  if (HaveNbrS)
  {
    flag = MPI_Wait(&reqRS, &stat);
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Wait = " << flag << endl;
      return -1;
    }
    flag = MPI_Wait(&reqSS, &stat);
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Wait = " << flag << endl;
      return -1;
    }
  }

  if (HaveNbrN)
  {
    flag = MPI_Wait(&reqRN, &stat);
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Wait = " << flag << endl;
      return -1;
    }
    flag = MPI_Wait(&reqSN, &stat);
    if (flag != MPI_SUCCESS)
    {
      cerr << "Error in MPI_Wait = " << flag << endl;
      return -1;
    }
  }

  // Return success
    return 0;
}


// -----------------------------------------------------------------------------
// UserData helper functions
// -----------------------------------------------------------------------------

UserData::~UserData()
{
  // Free exchange buffers
  free_buffers();

  // Free preconditioner data
  if (diag)
  {
    N_VDestroy(diag);
    diag = NULL;
  }
}


// -----------------------------------------------------------------------------
// UserOutput functions
// -----------------------------------------------------------------------------


// Parse command line inputs
int UserOutput::parse_args(vector<string> &args, bool outproc)
{
  vector<string>::iterator it;

  it = find(args.begin(), args.end(), "--help");
  if (it != args.end())
  {
    if (outproc) help();
    return 0;
  }

  it = find(args.begin(), args.end(), "--output");
  if (it != args.end())
  {
    output = stoi(*(it + 1));
    args.erase(it, it + 2);
  }

  it = find(args.begin(), args.end(), "--nout");
  if (it != args.end())
  {
    nout = stoi(*(it + 1));
    args.erase(it, it + 2);
  }

  return 0;
}


// Print command line options
void UserOutput::help()
{
  cout << endl;
  cout << "Output command line options:" << endl;
  cout << "  --output <level>  : output level" << endl;
  cout << "  --nout <nout>     : number of outputs" << endl;
}


void UserOutput::print()
{
  cout << endl;
  cout << " Output options:" << endl;
  cout << " --------------------------------- " << endl;
  cout << " output level = " << output << endl;
  cout << " nout         = " << nout   << endl;
  cout << " --------------------------------- " << endl;
}


int UserOutput::open(UserData *udata)
{
  bool outproc = (udata->myid_c == 0);

  // Header for status output
  if (outproc)
  {
    cout << scientific;
    cout << setprecision(numeric_limits<realtype>::digits10);
    cout << endl;
    if (error)
    {
      cout << "          t           ";
      cout << "          ||u||_rms      ";
      cout << "          max error      " << endl;
      cout << " ---------------------";
      cout << "-------------------------";
      cout << "-------------------------" << endl;
    }
    else
    {
      cout << "          t           ";
      cout << "          ||u||_rms      " << endl;
      cout << " ---------------------";
      cout << "-------------------------" << endl;
    }
  }

  // Output problem information and open output streams
  if (output == 2)
  {
    // Open output streams for solution and error
    stringstream fname;
    fname << "diffusion_2d_solution." << setfill('0') << setw(5)
          << udata->myid_c << ".txt";
    uoutstream.open(fname.str());

    uoutstream << "# title Diffusion 2D" << endl;
    uoutstream << "# nvar 1" << endl;
    uoutstream << "# vars u" << endl;
    uoutstream << "# nt  " << nout + 1   << endl;
    uoutstream << "# nx  " << udata->nx  << endl;
    uoutstream << "# xl  " << udata->xl  << endl;
    uoutstream << "# xu  " << udata->xu  << endl;
    uoutstream << "# ny  " << udata->ny  << endl;
    uoutstream << "# yl  " << udata->yl  << endl;
    uoutstream << "# yu  " << udata->yu  << endl;
    uoutstream << "# px  " << udata->npx << endl;
    uoutstream << "# py  " << udata->npy << endl;
    uoutstream << "# np  " << udata->np  << endl;
    uoutstream << "# is  " << udata->is  << endl;
    uoutstream << "# ie  " << udata->ie  << endl;
    uoutstream << "# js  " << udata->js  << endl;
    uoutstream << "# je  " << udata->je  << endl;

    uoutstream << scientific;
    uoutstream << setprecision(numeric_limits<realtype>::digits10);

    if (error)
    {
      fname.str("");
      fname.clear();
      fname << "diffusion_2d_error." << setfill('0') << setw(5) << udata->myid_c
            << ".txt";
      eoutstream.open(fname.str());

      eoutstream << "# title Diffusion 2D Error" << endl;
      eoutstream << "# nvar 1" << endl;
      eoutstream << "# vars u" << endl;
      eoutstream << "# nt  " << nout + 1   << endl;
      eoutstream << "# nx  " << udata->nx  << endl;
      eoutstream << "# xl  " << udata->xl  << endl;
      eoutstream << "# xu  " << udata->xu  << endl;
      eoutstream << "# ny  " << udata->ny  << endl;
      eoutstream << "# yl  " << udata->yl  << endl;
      eoutstream << "# yu  " << udata->yu  << endl;
      eoutstream << "# px  " << udata->npx << endl;
      eoutstream << "# py  " << udata->npy << endl;
      eoutstream << "# np  " << udata->np  << endl;
      eoutstream << "# is  " << udata->is  << endl;
      eoutstream << "# ie  " << udata->ie  << endl;
      eoutstream << "# js  " << udata->js  << endl;
      eoutstream << "# je  " << udata->je  << endl;

      eoutstream << scientific;
      eoutstream << setprecision(numeric_limits<realtype>::digits10);
    }
  }

  return 0;
}


int UserOutput::write(realtype t, N_Vector u, UserData *udata)
{
  int      flag;
  realtype max;
  bool     outproc = (udata->myid_c == 0);

  if (output > 0)
  {
    if (error)
    {
      // Compute the error
      flag = SolutionError(t, u, error, udata);
      if (check_flag(&flag, "SolutionError", 1)) return 1;

      // Compute max error
      max = N_VMaxNorm(error);
    }

    // Compute rms norm of the state
    realtype urms = sqrt(N_VDotProd(u, u) / udata->nx / udata->ny);

    // Output current status
    if (outproc)
    {
      if (error)
      {
        cout << setw(22) << t << setw(25) << urms << setw(25) << max << endl;
      }
      else
      {
        cout << setw(22) << t << setw(25) << urms << endl;
      }
    }

    // Write solution and error to disk
    if (output == 2)
    {
      // Sync host and device memory
      flag = CopyDataFromDevice(u);
      if (check_flag(&flag, "CopyDataFromDevice", 1)) return -1;

      realtype *uarray = N_VGetArrayPointer(u);
      if (check_flag((void *) uarray, "N_VGetArrayPointer", 0)) return -1;

      uoutstream << t << " ";
      for (sunindextype i = 0; i < udata->nodes_loc; i++)
      {
        uoutstream << uarray[i] << " ";
      }
      uoutstream << endl;

      if (error)
      {
        // Sync host and device memory
        flag = CopyDataFromDevice(error);
        if (check_flag(&flag, "CopyDataFromDevice", 1)) return -1;

        // Output error to disk
        realtype *earray = N_VGetArrayPointer(error);
        if (check_flag((void *) earray, "N_VGetArrayPointer", 0)) return -1;

        eoutstream << t << " ";
        for (sunindextype i = 0; i < udata->nodes_loc; i++)
        {
          eoutstream << earray[i] << " ";
        }
        eoutstream << endl;
      }
    }
  }
  return 0;
}


int UserOutput::close(UserData *udata)
{
  bool outproc = (udata->myid_c == 0);

  // Footer for status output
  if (outproc && (output > 0))
  {
    if (error)
    {
      cout << " ---------------------";
      cout << "-------------------------";
      cout << "-------------------------" << endl;
      cout << endl;
    }
    else
    {
      cout << " ---------------------";
      cout << "-------------------------" << endl;
      cout << endl;
    }
  }

  if (output == 2)
  {
    // Close output streams
    uoutstream.close();
    if (error) eoutstream.close();
  }

  if (error)
  {
    // Free error vector
    N_VDestroy(error);
    error = NULL;
  }

  return 0;
}


// -----------------------------------------------------------------------------
// Output and utility functions
// -----------------------------------------------------------------------------

// Compute the solution error
int SolutionError(realtype t, N_Vector u, N_Vector e, UserData *udata)
{
  // Compute true solution
  int flag = Solution(t, e, udata);
  if (flag != 0) return -1;

  // Compute absolute error
  N_VLinearSum(ONE, u, -ONE, e, e);
  N_VAbs(e, e);

  return 0;
}


// Check function return value
int check_flag(void *flagvalue, const string funcname, int opt)
{
  // Check if the function returned a NULL pointer
  if (opt == 0)
  {
    if (flagvalue == NULL)
    {
      cerr << endl << "ERROR: " << funcname << " returned NULL pointer" << endl
           << endl;
      return 1;
    }
  }
  // Check the function return flag value
  else if (opt == 1 || opt == 2)
  {
    int errflag = *((int *) flagvalue);
    if  ((opt == 1 && errflag < 0) || (opt == 2 && errflag != 0))
    {
      cerr << endl << "ERROR: " << funcname << " returned with flag = "
           << errflag << endl << endl;
      return 1;
    }
  }
  else
  {
    cerr << endl << "ERROR: check_flag called with an invalid option value"
         << endl;
    return 1;
  }

  return 0;
}

//---- end of file ----
