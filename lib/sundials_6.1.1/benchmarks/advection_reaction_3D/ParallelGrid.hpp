/* -----------------------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
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
 * A simple implementation of a parallel structured cartesian mesh class that
 * supports up to 3 dimensions and an arbitrary number of degrees of freedom.
 * ----------------------------------------------------------------------------*/

#ifndef _SIMPLEPARGRID_H
#define _SIMPLEPARGRID_H

#include <iomanip>
#include <iostream>
#include <fstream>
#include <mpi.h>

#include <sundials/sundials_memory.h>

namespace sundials_tools
{

// Types of boundaries supported.
enum class BoundaryType
{
  PERIODIC
};

// Types of stencils supported.
enum class StencilType
{
  UPWIND
};

template<typename REAL, typename GLOBALINT, int NDIMS>
class ParallelGrid
{
public:
  // Constructor that creates a new ParallelGrid object.
  // [in] - the memory helper to use for allocating the MPI buffers
  // [in,out] comm - on input, the overal MPI communicator, on output, the cartesian communicator
  // [in] a[] - an array of length NDIMS which defines the domain [a,b]
  // [in] b[] - an array of length NDIMS which defines the domain [a,b]
  // [in] npts[] - an array of length NDIMS which defines the number of mesh points in each dimension
  // [in] dof - the number of degrees of freedom in each dimension
  // [in] bc - the type of boundary conditions (see BoundaryType)
  // [in] st - the stencil to use (see StencilType)
  // [in] width - the stencil width; defaults to 1
  // [in] npxyz - the number of processors in each dimension; defaults to 0 which means MPI will choose
  // [in] reorder - should MPI_Cart_create do process reordering to optimize or not; defaults to false (some MPI implementations ignore this)
  ParallelGrid(SUNMemoryHelper memhelp, MPI_Comm* comm, const REAL a[], const REAL b[], const GLOBALINT npts[], int dof,
               BoundaryType bc, StencilType st, int width = 1, const int npxyz[] = nullptr, bool reorder = false)
    : nx(1), ny(1), nz(1),
      nxl(1), nyl(1), nzl(1),
      npx(1), npy(1), npz(1),
      dx(0.0), dy(0.0), dz(0.0),
      ax(0.0), ay(0.0), az(0.0),
      bx(0.0), by(0.0), bz(0.0),
      dof(dof), dims{0,0,0}, coords{0,0,0},
      bc(bc), st(st), width(width),
      memhelp(memhelp)
  {
    static_assert((NDIMS >= 1 && NDIMS <= 3), "ParallelGrid NDIMS must be 1, 2 or 3");

    int retval, nprocs;
    int periods[] = {0, 0, 0};

    if (npxyz)
    {
      dims[0] = npxyz[0];
      if (NDIMS >= 2) dims[1] = npxyz[1];
      if (NDIMS == 3) dims[2] = npxyz[2];
    }

    MPI_Comm_size(*comm, &nprocs);
    retval = MPI_Dims_create(nprocs, NDIMS, dims);
    assert(retval == MPI_SUCCESS);

    periods[0] = bc == BoundaryType::PERIODIC;
    periods[1] = bc == BoundaryType::PERIODIC;
    periods[2] = bc == BoundaryType::PERIODIC;
    retval = MPI_Cart_create(*comm, NDIMS, dims, periods, reorder, comm);
    assert(retval == MPI_SUCCESS);

    retval = MPI_Cart_get(*comm, NDIMS, dims, periods, coords);
    assert(retval == MPI_SUCCESS);

    cart_comm = *comm;

    npx = dims[0];
    nx  = npts[0];
    ax  = a[0];
    bx  = b[0];
    dx  = (bx-ax) / (REAL) nx;
    int is = nx*(coords[0])/npx;
    int ie = nx*(coords[0]+1)/npx-1;
    nxl = ie-is+1;

    neq = dof * nxl;

    if (NDIMS >= 2)
    {
      npy = dims[1];
      ny  = npts[1];
      ay  = a[1];
      by  = b[1];
      dy  = (by-ay) / (REAL) ny;
      int js = ny*(coords[1])/npy;
      int je = ny*(coords[1]+1)/npy-1;
      nyl = je-js+1;

      neq *= nyl;
    }

    if (NDIMS == 3)
    {
      npz = dims[2];
      nz  = npts[2];
      az  = a[2];
      bz  = b[2];
      dz  = (bz-az) / (REAL) nz;
      int ks = nz*(coords[2])/npz;
      int ke = nz*(coords[2]+1)/npz-1;
      nzl = ke-ks+1;

      neq *= nzl;
    }

    if (st == StencilType::UPWIND)
      AllocateBuffersUpwind();

  }

  // TODO:
  //  - does not take advantage of upwind scheme to reduce communications and memory
  //  - support non-periodic boundary conditions
  // For all faces where neighbors exist: determine neighbor process indices.
  // For all faces: allocate exchange buffers.
  void AllocateBuffersUpwind()
  {
    int retval = 0;
    int nbcoords[] = {0, 0, 0};

    SUNMemoryHelper_Alloc(memhelp, &Wrecv_, sizeof(REAL)*dof*width*nyl*nzl,
                          memoryType(), nullptr);
    SUNMemoryHelper_Alloc(memhelp, &Wsend_, sizeof(REAL)*dof*width*nyl*nzl,
                          memoryType(), nullptr);
    ipW = MPI_PROC_NULL;
    if ((coords[0] > 0) || (bc == BoundaryType::PERIODIC)) {
      nbcoords[0] = coords[0]-1;
      nbcoords[1] = coords[1];
      nbcoords[2] = coords[2];
      retval = MPI_Cart_rank(cart_comm, nbcoords, &ipW);
      assert(retval == MPI_SUCCESS);
    }

    SUNMemoryHelper_Alloc(memhelp, &Erecv_, sizeof(REAL)*dof*width*nyl*nzl,
                          memoryType(), nullptr);
    SUNMemoryHelper_Alloc(memhelp, &Esend_, sizeof(REAL)*dof*width*nyl*nzl,
                          memoryType(), nullptr);
    ipE = MPI_PROC_NULL;
    if ((coords[0] < dims[0]-1) || (bc == BoundaryType::PERIODIC)) {
      nbcoords[0] = coords[0]+1;
      nbcoords[1] = coords[1];
      nbcoords[2] = coords[2];
      retval = MPI_Cart_rank(cart_comm, nbcoords, &ipE);
      assert(retval == MPI_SUCCESS);
    }

    if (NDIMS >= 2)
    {
      SUNMemoryHelper_Alloc(memhelp, &Srecv_, sizeof(REAL)*dof*width*nxl*nzl,
                            memoryType(), nullptr);
      SUNMemoryHelper_Alloc(memhelp, &Ssend_, sizeof(REAL)*dof*width*nxl*nzl,
                            memoryType(), nullptr);
      ipS = MPI_PROC_NULL;
      if ((coords[1] > 0) || (bc == BoundaryType::PERIODIC)) {
        nbcoords[0] = coords[0];
        nbcoords[1] = coords[1]-1;
        nbcoords[2] = coords[2];
        retval = MPI_Cart_rank(cart_comm, nbcoords, &ipS);
        assert(retval == MPI_SUCCESS);
      }

      SUNMemoryHelper_Alloc(memhelp, &Nrecv_, sizeof(REAL)*dof*width*nxl*nzl,
                            memoryType(), nullptr);
      SUNMemoryHelper_Alloc(memhelp, &Nsend_, sizeof(REAL)*dof*width*nxl*nzl,
                            memoryType(), nullptr);
      ipN = MPI_PROC_NULL;
      if ((coords[1] < dims[1]-1) || (bc == BoundaryType::PERIODIC)) {
        nbcoords[0] = coords[0];
        nbcoords[1] = coords[1]+1;
        nbcoords[2] = coords[2];
        retval = MPI_Cart_rank(cart_comm, nbcoords, &ipN);
        assert(retval == MPI_SUCCESS);
      }
    }

    if (NDIMS == 3)
    {
      SUNMemoryHelper_Alloc(memhelp, &Brecv_, sizeof(REAL)*dof*width*nxl*nyl,
                            memoryType(), nullptr);
      SUNMemoryHelper_Alloc(memhelp, &Bsend_, sizeof(REAL)*dof*width*nxl*nyl,
                            memoryType(), nullptr);
      ipB = MPI_PROC_NULL;
      if ((coords[2] > 0) || (bc == BoundaryType::PERIODIC)) {
        nbcoords[0] = coords[0];
        nbcoords[1] = coords[1];
        nbcoords[2] = coords[2]-1;
        retval = MPI_Cart_rank(cart_comm, nbcoords, &ipB);
        assert(retval == MPI_SUCCESS);
      }

      SUNMemoryHelper_Alloc(memhelp, &Frecv_, sizeof(REAL)*dof*width*nxl*nyl,
                            memoryType(), nullptr);
      SUNMemoryHelper_Alloc(memhelp, &Fsend_, sizeof(REAL)*dof*width*nxl*nyl,
                            memoryType(), nullptr);
      ipF = MPI_PROC_NULL;
      if ((coords[2] < dims[2]-1) || (bc == BoundaryType::PERIODIC)) {
        nbcoords[0] = coords[0];
        nbcoords[1] = coords[1];
        nbcoords[2] = coords[2]+1;
        retval = MPI_Cart_rank(cart_comm, nbcoords, &ipF);
        assert(retval == MPI_SUCCESS);
      }
    }

  }

  // TODO: this could be optimized for upwind
  int ExchangeStart(std::function<void (REAL*,REAL*,REAL*,REAL*,REAL*,REAL*)> fill)
  {
    int retval = 0;

    // Initialize all requests in array
    for (int i=0; i<12; i++)
      req[i] = MPI_REQUEST_NULL;

    // Open an Irecv buffer for each neighbor
    if (ipW != MPI_PROC_NULL)
    {
      retval = MPI_Irecv(getRecvBuffer("EAST"), dof*nyl*nzl, MPI_SUNREALTYPE, ipW,
                         1, cart_comm, req);
      assert(retval == MPI_SUCCESS);
    }

    if (ipE != MPI_PROC_NULL)
    {
      retval = MPI_Irecv(getRecvBuffer("WEST"), dof*nyl*nzl, MPI_SUNREALTYPE, ipE,
                         0, cart_comm, req+1);
      assert(retval == MPI_SUCCESS);
    }

    if (NDIMS >= 2)
    {
      if (ipS != MPI_PROC_NULL)
      {
        retval = MPI_Irecv(getRecvBuffer("NORTH"), dof*nxl*nzl, MPI_SUNREALTYPE, ipS,
                           3, cart_comm, req+2);
        assert(retval == MPI_SUCCESS);
      }

      if (ipN != MPI_PROC_NULL)
      {
        retval = MPI_Irecv(getRecvBuffer("SOUTH"), dof*nxl*nzl, MPI_SUNREALTYPE, ipN,
                           2, cart_comm, req+3);
        assert(retval == MPI_SUCCESS);
      }
    }

    if (NDIMS >= 3)
    {
      if (ipB != MPI_PROC_NULL)
      {
        retval = MPI_Irecv(getRecvBuffer("FRONT"), dof*nxl*nyl, MPI_SUNREALTYPE, ipB,
                           5, cart_comm, req+4);
        assert(retval == MPI_SUCCESS);
      }

      if (ipF != MPI_PROC_NULL)
      {
        retval = MPI_Irecv(getRecvBuffer("BACK"), dof*nxl*nyl, MPI_SUNREALTYPE, ipF,
                           4, cart_comm, req+5);
        assert(retval == MPI_SUCCESS);
      }
    }

    // Call user lambda to fill the send buffers
    fill(getSendBuffer("WEST"),
         getSendBuffer("EAST"),
         getSendBuffer("SOUTH"),
         getSendBuffer("NORTH"),
         getSendBuffer("BACK"),
         getSendBuffer("FRONT"));

    // Send data to neighbors
    if (ipW != MPI_PROC_NULL)
    {
      retval = MPI_Isend(getSendBuffer("EAST"), dof*nyl*nzl, MPI_SUNREALTYPE, ipW, 0,
                         cart_comm, req+6);
      assert(retval == MPI_SUCCESS);
    }

    if (ipE != MPI_PROC_NULL)
    {
      retval = MPI_Isend(getSendBuffer("WEST"), dof*nyl*nzl, MPI_SUNREALTYPE, ipE, 1,
                         cart_comm, req+7);
      assert(retval == MPI_SUCCESS);
    }

    if (NDIMS >= 2)
    {
      if (ipS != MPI_PROC_NULL)
      {
        retval = MPI_Isend(getSendBuffer("NORTH"), dof*nxl*nzl, MPI_SUNREALTYPE, ipS, 2,
                           cart_comm, req+8);
        assert(retval == MPI_SUCCESS);
      }

      if (ipN != MPI_PROC_NULL)
      {
        retval = MPI_Isend(getSendBuffer("SOUTH"), dof*nxl*nzl, MPI_SUNREALTYPE, ipN, 3,
                           cart_comm, req+9);
        assert(retval == MPI_SUCCESS);
      }
    }

    if (NDIMS == 3)
    {
      if (ipB != MPI_PROC_NULL)
      {
        retval = MPI_Isend(getSendBuffer("FRONT"), dof*nxl*nyl, MPI_SUNREALTYPE, ipB, 4,
                           cart_comm, req+10);
        assert(retval == MPI_SUCCESS);
      }

      if (ipF != MPI_PROC_NULL)
      {
        retval = MPI_Isend(getSendBuffer("BACK"), dof*nxl*nyl, MPI_SUNREALTYPE, ipF, 5,
                           cart_comm, req+11);
        assert(retval == MPI_SUCCESS);
      }
    }

    return retval;
  }

  // Waits for neighbor exchange to finish.
  int ExchangeEnd()
  {
    MPI_Status stat[12];
    int retval;

    // Wait for messages to finish send/receive
    retval = MPI_Waitall(12, req, stat);
    assert(retval == MPI_SUCCESS);

    return retval;
  }

  // Prints out information about the ParallelGrid to stdout.
  void PrintInfo()
  {
    printf("ParallelGrid Info:\n");
    printf("    dimensions = %d\n", NDIMS);
    printf("    processors = {%d, %d, %d}\n", npx, npy, npz);
    printf("        domain = {[%g,%g], [%g,%g], [%g,%g]}\n", ax, bx, ay, by, az, bz);
    printf("   global npts = {%li, %li, %li}\n", (long int) nx, (long int) ny, (long int) nz);
    printf("    local npts = {%d, %d, %d}\n", nxl, nyl, nzl);
    printf("  mesh spacing = {%g, %g, %g}\n", dx, dy, dz);
  }

  // Saves the mesh to a file.
  //    First row is x. Second row is y. Third row is z.
  //    Can be loaded into MATLAB like so:
  //      mesh = loadtxt('mesh.txt');
  //      [X,Y,Z] = meshgrid(mesh(1,:),mesh(2,:),mesh(3,:));
  void MeshToFile(const std::string& fname)
  {
    std::ofstream mesh_file;
    mesh_file.open(fname);
    mesh_file << std::setprecision(16);
    for (GLOBALINT i = 0; i < nx; i++)
      mesh_file << " " << dx*i;
    mesh_file << std::endl;
    for (GLOBALINT i = 0; i < ny; i++)
      mesh_file << " " << dy*i;
    mesh_file << std::endl;
    for (GLOBALINT i = 0; i < nz; i++)
      mesh_file << " " << dz*i;
    mesh_file << std::endl;
    mesh_file.close();
  }

  int nprocs() const
  {
    return npx*npy*npz;
  }

  GLOBALINT npts() const
  {
    if (NDIMS == 1) return nx;
    if (NDIMS == 2) return nx*ny;
    if (NDIMS == 3) return nx*ny*nz;
  }

  GLOBALINT nptsl() const
  {
    if (NDIMS == 1) return nxl;
    if (NDIMS == 2) return nxl*nyl;
    if (NDIMS == 3) return nxl*nyl*nzl;
  }

  GLOBALINT neql() const
  {
    return dof*nptsl();
  }

  REAL* getRecvBuffer(const std::string& direction)
  {
    if (direction == "WEST")
    {
      return static_cast<REAL*>(Wrecv_->ptr);
    }
    else if (direction == "EAST")
    {
      return static_cast<REAL*>(Erecv_->ptr);
    }
    else if (direction == "NORTH")
    {
      return static_cast<REAL*>(Nrecv_->ptr);
    }
    else if (direction == "SOUTH")
    {
      return static_cast<REAL*>(Srecv_->ptr);
    }
    else if (direction == "FRONT")
    {
      return static_cast<REAL*>(Frecv_->ptr);
    }
    else if (direction == "BACK")
    {
      return static_cast<REAL*>(Brecv_->ptr);
    }
    else
    {
      return nullptr;
    }
  }

  REAL* getSendBuffer(const std::string& direction)
  {
    if (direction == "WEST")
    {
      return static_cast<REAL*>(Wsend_->ptr);
    }
    else if (direction == "EAST")
    {
      return static_cast<REAL*>(Esend_->ptr);
    }
    else if (direction == "NORTH")
    {
      return static_cast<REAL*>(Nsend_->ptr);
    }
    else if (direction == "SOUTH")
    {
      return static_cast<REAL*>(Ssend_->ptr);
    }
    else if (direction == "FRONT")
    {
      return static_cast<REAL*>(Fsend_->ptr);
    }
    else if (direction == "BACK")
    {
      return static_cast<REAL*>(Bsend_->ptr);
    }
    else
    {
      return nullptr;
    }
  }

  ~ParallelGrid()
  {
    SUNMemoryHelper_Dealloc(memhelp, Esend_, nullptr);
    SUNMemoryHelper_Dealloc(memhelp, Wsend_, nullptr);
    SUNMemoryHelper_Dealloc(memhelp, Nsend_, nullptr);
    SUNMemoryHelper_Dealloc(memhelp, Ssend_, nullptr);
    SUNMemoryHelper_Dealloc(memhelp, Fsend_, nullptr);
    SUNMemoryHelper_Dealloc(memhelp, Bsend_, nullptr);
    SUNMemoryHelper_Dealloc(memhelp, Erecv_, nullptr);
    SUNMemoryHelper_Dealloc(memhelp, Wrecv_, nullptr);
    SUNMemoryHelper_Dealloc(memhelp, Nrecv_, nullptr);
    SUNMemoryHelper_Dealloc(memhelp, Srecv_, nullptr);
    SUNMemoryHelper_Dealloc(memhelp, Frecv_, nullptr);
    SUNMemoryHelper_Dealloc(memhelp, Brecv_, nullptr);
  }

  GLOBALINT nx, ny, nz;    /* number of intervals globally       */
  int       nxl, nyl, nzl; /* number of intervals locally        */
  int       npx, npy, npz; /* numner of processes                */
  REAL      dx, dy, dz;    /* mesh spacing                       */
  REAL      ax, ay, az;    /* domain in [a, b]                   */
  REAL      bx, by, bz;
  int       dof;           /* degrees of freedom per node        */
  int       neq;           /* total number of equations locally  */

  int       ipW, ipE;      /* MPI ranks for neighbor procs       */
  int       ipS, ipN;
  int       ipB, ipF;

  int       dims[3];
  int       coords[3];


private:
  MPI_Comm     cart_comm;  /* MPI cartesian communicator         */
  MPI_Request  req[12];

  BoundaryType bc;
  StencilType  st;
  int          width;

  SUNMemoryHelper memhelp;
  SUNMemory Wsend_;            /* MPI send/recv buffers              */
  SUNMemory Esend_;
  SUNMemory Ssend_;
  SUNMemory Nsend_;
  SUNMemory Bsend_;
  SUNMemory Fsend_;
  SUNMemory Wrecv_;
  SUNMemory Erecv_;
  SUNMemory Srecv_;
  SUNMemory Nrecv_;
  SUNMemory Brecv_;
  SUNMemory Frecv_;

  SUNMemoryType memoryType()
  {
    SUNMemory test;
    if (!SUNMemoryHelper_Alloc(memhelp, &test, sizeof(REAL), SUNMEMTYPE_PINNED,
                               nullptr))
    {
      SUNMemoryHelper_Dealloc(memhelp, test, nullptr);
      return(SUNMEMTYPE_PINNED);
    }
    if (!SUNMemoryHelper_Alloc(memhelp, &test, sizeof(REAL), SUNMEMTYPE_DEVICE,
                               nullptr))
    {
      SUNMemoryHelper_Dealloc(memhelp, test, nullptr);
      return(SUNMEMTYPE_DEVICE);
    }
    if (!SUNMemoryHelper_Alloc(memhelp, &test, sizeof(REAL), SUNMEMTYPE_UVM,
                               nullptr))
    {
      SUNMemoryHelper_Dealloc(memhelp, test, nullptr);
      return(SUNMEMTYPE_UVM);
    }
    else
    {
      return(SUNMEMTYPE_HOST);
    }
  }

};

}

#endif
