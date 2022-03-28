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
 * Serial functions for exchange buffers
 * ---------------------------------------------------------------------------*/

#include "diffusion_2D.hpp"

// Pack exchange buffers
int UserData::pack_buffers(const N_Vector u)
{
  // Access data array
  const realtype *uarray = N_VGetArrayPointer(u);
  if (check_flag((void *) uarray, "N_VGetArrayPointer", 0)) return -1;

  if (HaveNbrW)
    for (sunindextype i = 0; i < ny_loc; i++)
      Wsend[i] = uarray[IDX(0,i,nx_loc)];

  if (HaveNbrE)
    for (sunindextype i = 0; i < ny_loc; i++)
      Esend[i] = uarray[IDX(nx_loc-1,i,nx_loc)];

  if (HaveNbrS)
    for (sunindextype i = 0; i < nx_loc; i++)
      Ssend[i] = uarray[IDX(i,0,nx_loc)];

  if (HaveNbrN)
    for (sunindextype i = 0; i < nx_loc; i++)
      Nsend[i] = uarray[IDX(i,ny_loc-1,nx_loc)];

  return 0;
}


// Allocate exchange buffers
int UserData::allocate_buffers()
{
  if (HaveNbrW)
  {
    Wrecv = new realtype[ny_loc];
    Wsend = new realtype[ny_loc];
  }

  if (HaveNbrE)
  {
    Erecv = new realtype[ny_loc];
    Esend = new realtype[ny_loc];
  }

  if (HaveNbrS)
  {
    Srecv = new realtype[nx_loc];
    Ssend = new realtype[nx_loc];
  }

  if (HaveNbrN)
  {
    Nrecv = new realtype[nx_loc];
    Nsend = new realtype[nx_loc];
  }

  return 0;
}


// Free exchange buffers
int UserData::free_buffers()
{
  // Free exchange buffers
  if (Wrecv != NULL)  delete[] Wrecv;
  if (Wsend != NULL)  delete[] Wsend;
  if (Erecv != NULL)  delete[] Erecv;
  if (Esend != NULL)  delete[] Esend;
  if (Srecv != NULL)  delete[] Srecv;
  if (Ssend != NULL)  delete[] Ssend;
  if (Nrecv != NULL)  delete[] Nrecv;
  if (Nsend != NULL)  delete[] Nsend;

  return 0;
}
