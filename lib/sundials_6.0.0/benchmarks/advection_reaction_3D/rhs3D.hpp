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
 * -----------------------------------------------------------------------------*/

#ifndef ADVECTION_REACTION_3D_RHS_HPP
#define ADVECTION_REACTION_3D_RHS_HPP

#include "advection_reaction_3D.hpp"

using raja_xyz_tuple = camp::tuple<RAJA::RangeSegment, RAJA::RangeSegment, RAJA::RangeSegment>;

/* --------------------------------------------------------------
 * Right hand side (RHS) and residual functions
 * --------------------------------------------------------------*/

/* Compute the advection term f(t,y) = -c (grad * y). This is done using
   upwind 1st order finite differences. */
static int Advection(realtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  /* access problem data */
  UserData* udata = (UserData*) user_data;

  SUNDIALS_CXX_MARK_FUNCTION(udata->prof);

  /* set variable shortcuts */
  const int      nxl = udata->grid->nxl;
  const int      nyl = udata->grid->nyl;
  const int      nzl = udata->grid->nzl;
  const int      dof = udata->grid->dof;
  const realtype c   = udata->c;
  const realtype cx  = -c / udata->grid->dx;
  const realtype cy  = -c / udata->grid->dy;
  const realtype cz  = -c / udata->grid->dz;

  /* local variables */
  int retval;

  /* begin exchanging boundary information */
  if (udata->grid->nprocs() > 1)
  {
    retval = ExchangeAllStart(y, udata);
    if (check_retval(&retval, "ExchangeAllStart", 1, udata->myid))
      return(-1);
  }

  /* set output to zero */
  N_VConst(0.0, ydot);

  /* create views of the data */
  RAJA::View<realtype, RAJA::Layout<NDIMS+1> > Yview(GetVecData(y),
                                                     nxl, nyl, nzl, dof);
  RAJA::View<realtype, RAJA::Layout<NDIMS+1> > dYview(GetVecData(ydot),
                                                      nxl, nyl, nzl, dof);

  /* iterate over domain interior, computing advection */
  if (c > 0.0)
  {
    /* flow moving in the positive x,y,z direction */
    auto range = RAJA::make_tuple(RAJA::RangeSegment(1, nxl),
                                  RAJA::RangeSegment(1, nyl),
                                  RAJA::RangeSegment(1, nzl));

    RAJA::kernel<XYZ_KERNEL_POL>(range,
      [=] DEVICE_FUNC (int i, int j, int k) {
      const realtype u_ijk = Yview(i,j,k,0);
      const realtype v_ijk = Yview(i,j,k,1);
      const realtype w_ijk = Yview(i,j,k,2);

      // grad * u
      dYview(i,j,k,0) =  cz * (u_ijk - Yview(i,j,k-1,0)); // du/dz
      dYview(i,j,k,0) += cy * (u_ijk - Yview(i,j-1,k,0)); // du/dy
      dYview(i,j,k,0) += cx * (u_ijk - Yview(i-1,j,k,0)); // du/dx

      // grad * v
      dYview(i,j,k,1) =  cz * (v_ijk - Yview(i,j,k-1,1)); // dv/dz
      dYview(i,j,k,1) += cy * (v_ijk - Yview(i,j-1,k,1)); // dv/dy
      dYview(i,j,k,1) += cx * (v_ijk - Yview(i-1,j,k,1)); // dv/dx

      // grad * w
      dYview(i,j,k,2) =  cz * (w_ijk - Yview(i,j,k-1,2)); // dw/dz
      dYview(i,j,k,2) += cy * (w_ijk - Yview(i,j-1,k,2)); // dw/dy
      dYview(i,j,k,2) += cx * (w_ijk - Yview(i-1,j,k,2)); // dw/dx
    });
  }
  else if (c < 0.0)
  {
    /* flow moving in the negative x,y,z direction */
    auto range = RAJA::make_tuple(RAJA::RangeSegment(0, nxl-1),
                                  RAJA::RangeSegment(0, nyl-1),
                                  RAJA::RangeSegment(0, nzl-1));
    RAJA::kernel<XYZ_KERNEL_POL>(range,
      [=] DEVICE_FUNC (int i, int j, int k) {
      const realtype u_ijk = Yview(i,j,k,0);
      const realtype v_ijk = Yview(i,j,k,1);
      const realtype w_ijk = Yview(i,j,k,2);

      // grad * u
      dYview(i,j,k,0) =  cz * (u_ijk - Yview(i,j,k+1,0)); // du/dz
      dYview(i,j,k,0) += cy * (u_ijk - Yview(i,j+1,k,0)); // du/dy
      dYview(i,j,k,0) += cx * (u_ijk - Yview(i+1,j,k,0)); // du/dx

      // grad * v
      dYview(i,j,k,1) =  cz * (v_ijk - Yview(i,j,k+1,1)); // dv/dz
      dYview(i,j,k,1) += cy * (v_ijk - Yview(i,j+1,k,1)); // dv/dy
      dYview(i,j,k,1) += cx * (v_ijk - Yview(i+1,j,k,1)); // dv/dx

      // grad * w
      dYview(i,j,k,2) =  cz * (w_ijk - Yview(i,j,k+1,2)); // dw/dz
      dYview(i,j,k,2) += cy * (w_ijk - Yview(i,j+1,k,2)); // dw/dy
      dYview(i,j,k,2) += cx * (w_ijk - Yview(i+1,j,k,2)); // dw/dx
    });
  }

  /* finish exchanging boundary information */
  if (udata->grid->nprocs() > 1)
  {
    retval = ExchangeAllEnd(udata);
    if (check_retval(&retval, "ExchangeAllEnd", 1, udata->myid))
      return(-1);
  }

  /* compute advection at process boundaries */
  if (c > 0.0)
  {
    if (udata->grid->npx > 1)
    {
      /* Flow moving in the positive x,y,z direction:
      *  boundaries are west face, south face, front face */

      RAJA::View<realtype, RAJA::Layout<NDIMS> >
        Yim1jk(udata->grid->getRecvBuffer("WEST"), nyl, nzl, dof); // Wrecv should have data that was sent from East

      auto west_face = RAJA::make_tuple(RAJA::RangeSegment(0, nyl),
                                        RAJA::RangeSegment(0, nzl),
                                        RAJA::RangeSegment(0, dof));

      RAJA::kernel<XYZ_KERNEL_POL>(west_face,
        [=] DEVICE_FUNC (int j, int k, int l) {
        dYview(0,j,k,l) += cx * (Yview(0,j,k,l) - Yim1jk(j,k,l)); // d/dx
      });
    }
    else
    {
      auto range = RAJA::make_tuple(RAJA::RangeSegment(0, 1),
                                    RAJA::RangeSegment(0, 1),
                                    RAJA::RangeSegment(0, 1));

      RAJA::kernel<XYZ_KERNEL_POL>(range,
        [=] DEVICE_FUNC (int i, int j, int k) {
        const realtype u_ijk = Yview(i,j,k,0);
        const realtype v_ijk = Yview(i,j,k,1);
        const realtype w_ijk = Yview(i,j,k,2);

        dYview(i,j,k,0) = cx * (u_ijk - Yview(nxl-1,j,k,0)); // du/dx
        dYview(i,j,k,1) = cx * (v_ijk - Yview(nxl-1,j,k,1)); // dv/dx
        dYview(i,j,k,2) = cx * (w_ijk - Yview(nxl-1,j,k,2)); // dw/dx
      });

    }

    if (udata->grid->npy > 1)
    {
      RAJA::View<realtype, RAJA::Layout<NDIMS> >
        Yijm1k(udata->grid->getRecvBuffer("SOUTH"), nxl, nzl, dof); // Nrecv should have data that was sent from North

      auto south_face = RAJA::make_tuple(RAJA::RangeSegment(0, nxl),
                                         RAJA::RangeSegment(0, nzl),
                                         RAJA::RangeSegment(0, dof));

      RAJA::kernel<XYZ_KERNEL_POL>(south_face,
        [=] DEVICE_FUNC (int i, int k, int l) {
        dYview(i,0,k,l) += cy * (Yview(i,0,k,l) - Yijm1k(i,k,l)); // d/dy
      });
    }
    else
    {
      auto range = RAJA::make_tuple(RAJA::RangeSegment(0, 1),
                                    RAJA::RangeSegment(0, 1),
                                    RAJA::RangeSegment(0, 1));

      RAJA::kernel<XYZ_KERNEL_POL>(range,
        [=] DEVICE_FUNC (int i, int j, int k) {
        const realtype u_ijk = Yview(i,j,k,0);
        const realtype v_ijk = Yview(i,j,k,1);
        const realtype w_ijk = Yview(i,j,k,2);

        dYview(i,j,k,0) += cy * (u_ijk - Yview(i,nyl-1,k,0)); // du/dy
        dYview(i,j,k,1) += cy * (v_ijk - Yview(i,nyl-1,k,1)); // dv/dy
        dYview(i,j,k,2) += cy * (w_ijk - Yview(i,nyl-1,k,2)); // dw/dy
      });
    }

    if (udata->grid->npz > 1)
    {
      RAJA::View<realtype, RAJA::Layout<NDIMS> >
        Yijkm1(udata->grid->getRecvBuffer("FRONT"), nxl, nyl, dof); // Frecv should have data that was sent from Back

      auto front_face = RAJA::make_tuple(RAJA::RangeSegment(0, nxl),
                                         RAJA::RangeSegment(0, nyl),
                                         RAJA::RangeSegment(0, dof));

      RAJA::kernel<XYZ_KERNEL_POL>(front_face,
        [=] DEVICE_FUNC (int i, int j, int l) {
        dYview(i,j,0,l) += cz * (Yview(i,j,0,l) - Yijkm1(i,j,l)); // d/dz
      });

    }
    else
    {
      auto range = RAJA::make_tuple(RAJA::RangeSegment(0, 1),
                                    RAJA::RangeSegment(0, 1),
                                    RAJA::RangeSegment(0, 1));

      RAJA::kernel<XYZ_KERNEL_POL>(range,
        [=] DEVICE_FUNC (int i, int j, int k) {
        const realtype u_ijk = Yview(i,j,k,0);
        const realtype v_ijk = Yview(i,j,k,1);
        const realtype w_ijk = Yview(i,j,k,2);

        dYview(i,j,k,0) +=  cz * (u_ijk - Yview(i,j,nzl-1,0)); // du/dz
        dYview(i,j,k,1) +=  cz * (v_ijk - Yview(i,j,nzl-1,1)); // dv/dz
        dYview(i,j,k,2) +=  cz * (w_ijk - Yview(i,j,nzl-1,2)); // dw/dz
      });
    }
  }
  else if (c < 0.0)
  {
    if (udata->grid->nprocs() != 1)
    {
      /* Flow moving in the negative x,y,z direction:
      *  boundaries are west face, south face, and front face */

      RAJA::View<realtype, RAJA::Layout<3> >
        Yip1jk(udata->grid->getRecvBuffer("EAST"), nyl, nzl, dof);
      RAJA::View<realtype, RAJA::Layout<3> >
        Yijp1k(udata->grid->getRecvBuffer("NORTH"), nxl, nzl, dof);
      RAJA::View<realtype, RAJA::Layout<3> >
        Yijkp1(udata->grid->getRecvBuffer("BACK"), nxl, nyl, dof);

      auto front_face = RAJA::make_tuple(RAJA::RangeSegment(0, nxl-1),
                                         RAJA::RangeSegment(0, nyl-1),
                                         RAJA::RangeSegment(0, dof));
      RAJA::kernel<XYZ_KERNEL_POL>(front_face,
        [=] DEVICE_FUNC (int i, int j, int l) {
        dYview(i,j,0,l) =  cz * (Yview(i,j,0,l) - Yijkp1(i,nzl+1,l)); // d/dz
        dYview(i,j,0,l) += cy * (Yview(i,j,0,l) - Yijp1k(0,j+1,l));   // d/dy
        dYview(i,j,0,l) += cx * (Yview(i,j,0,l) - Yip1jk(i+1,0,l));   // d/dx
      });

      auto south_face = RAJA::make_tuple(RAJA::RangeSegment(0, nxl-1),
                                         RAJA::RangeSegment(0, nzl-1),
                                         RAJA::RangeSegment(0, dof));
      RAJA::kernel<XYZ_KERNEL_POL>(south_face,
        [=] DEVICE_FUNC (int i, int k, int l) {
        dYview(i,0,k,l) =  cz * (Yview(i,0,k,l) - Yijkp1(i,k+1,l));   // d/dz
        dYview(i,0,k,l) += cy * (Yview(i,0,k,l) - Yijp1k(0,nyl+1,l)); // d/dy
        dYview(i,0,k,l) += cx * (Yview(i,0,k,l) - Yip1jk(i+1,0,l));   // d/dx
      });

      auto east_face = RAJA::make_tuple(RAJA::RangeSegment(0, nyl-1),
                                        RAJA::RangeSegment(0, nzl-1),
                                        RAJA::RangeSegment(0, dof));
      RAJA::kernel<XYZ_KERNEL_POL>(east_face,
        [=] DEVICE_FUNC (int j, int k, int l) {
        dYview(0,j,k,l) =  cz * (Yview(0,j,k,l) - Yijkp1(0,k+1,l));   // d/dz
        dYview(0,j,k,l) += cy * (Yview(0,j,k,l) - Yijp1k(0,j+1,l));   // d/dy
        dYview(0,j,k,l) += cx * (Yview(0,j,k,l) - Yip1jk(nxl+1,0,l)); // d/dx
      });
    }
    else
    {
      auto range = RAJA::make_tuple(RAJA::RangeSegment(nxl-2, nxl),
                                    RAJA::RangeSegment(nyl-2, nyl),
                                    RAJA::RangeSegment(nzl-2, nzl));
      RAJA::kernel<XYZ_KERNEL_POL>(range,
        [=] DEVICE_FUNC (int i, int j, int k) {
        const realtype u_ijk = Yview(i,j,k,0);
        const realtype v_ijk = Yview(i,j,k,1);
        const realtype w_ijk = Yview(i,j,k,2);

        // grad * u
        dYview(i,j,k,0) =  cz * (u_ijk - Yview(i,j,0,0)); // du/dz
        dYview(i,j,k,0) += cy * (u_ijk - Yview(i,0,k,0)); // du/dy
        dYview(i,j,k,0) += cx * (u_ijk - Yview(0,j,k,0)); // du/dx

        // grad * v
        dYview(i,j,k,1) =  cz * (v_ijk - Yview(i,j,0,1)); // dv/dz
        dYview(i,j,k,1) += cy * (v_ijk - Yview(i,0,k,1)); // dv/dy
        dYview(i,j,k,1) += cx * (v_ijk - Yview(0,j,k,1)); // dv/dx

        // grad * w
        dYview(i,j,k,2) =  cz * (w_ijk - Yview(i,j,0,2)); // dw/dz
        dYview(i,j,k,2) += cy * (w_ijk - Yview(i,0,k,2)); // dw/dy
        dYview(i,j,k,2) += cx * (w_ijk - Yview(0,j,k,2)); // dw/dx
      });
    }
  }

  /* return success */
  return(0);
}


/* Compute the reaction term g(t,y). */
static int Reaction(realtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  /* access problem data */
  UserData* udata = (UserData*) user_data;

  SUNDIALS_CXX_MARK_FUNCTION(udata->prof);

  /* set variable shortcuts */
  const realtype A  = udata->A;
  const realtype B  = udata->B;
  const realtype k1 = udata->k1;
  const realtype k2 = udata->k2;
  const realtype k3 = udata->k3;
  const realtype k4 = udata->k4;
  const realtype k5 = udata->k5;
  const realtype k6 = udata->k6;

  /* local variables */
  realtype* Ydata  = NULL;
  realtype* dYdata = NULL;

  /* access data arrays */
  Ydata = GetVecData(y);
  if (check_retval((void *)Ydata, "GetVecData", 0, udata->myid))
    return(-1);

  dYdata = GetVecData(ydot);
  if (check_retval((void *)dYdata, "GetVecData", 0, udata->myid))
    return(-1);

  RAJA::View<realtype, RAJA::Layout<NDIMS+1> > Yview(GetVecData(y),
                                                     udata->grid->nxl,
                                                     udata->grid->nyl,
                                                     udata->grid->nzl,
                                                     udata->grid->dof);

  RAJA::View<realtype, RAJA::Layout<NDIMS+1> > dYview(GetVecData(ydot),
                                                      udata->grid->nxl,
                                                      udata->grid->nyl,
                                                      udata->grid->nzl,
                                                      udata->grid->dof);

  auto range = RAJA::make_tuple(RAJA::RangeSegment(0, udata->grid->nxl),
                                RAJA::RangeSegment(0, udata->grid->nyl),
                                RAJA::RangeSegment(0, udata->grid->nzl));

  /* iterate over domain, computing reactions */
  if (udata->add_reactions)
  {
    /* when we are not additively splitting the rhs, we add to ydot
       as we expect it to hold the advection term already */
    RAJA::kernel<XYZ_KERNEL_POL>(range,
      [=] DEVICE_FUNC (int i, int j, int k) {
      const realtype u = Yview(i,j,k,0);
      const realtype v = Yview(i,j,k,1);
      const realtype w = Yview(i,j,k,2);
      dYview(i,j,k,0) += k1 * A - k2 * w * u + k3 * u * u * v - k4 * u;
      dYview(i,j,k,1) += k2 * w * u - k3 * u * u * v;
      dYview(i,j,k,2) += -k2 * w * u + k5 * B - k6 * w;
    });
  }
  else
  {
    /* set output to zero */
    N_VConst(0.0, ydot);

    RAJA::kernel<XYZ_KERNEL_POL>(range,
      [=] DEVICE_FUNC (int i, int j, int k) {
      const realtype u = Yview(i,j,k,0);
      const realtype v = Yview(i,j,k,1);
      const realtype w = Yview(i,j,k,2);
      dYview(i,j,k,0) = k1 * A - k2 * w * u + k3 * u * u * v - k4 * u;
      dYview(i,j,k,1) = k2 * w * u - k3 * u * u * v;
      dYview(i,j,k,2) = -k2 * w * u + k5 * B - k6 * w;
    });
  }

  /* return success */
    return(0);
}


/* Compute the RHS as h(t,y) = f(t,y) + g(t,y). */
static int AdvectionReaction(realtype t, N_Vector y, N_Vector ydot,
                             void *user_data)
{
  /* access problem data */
  UserData* udata = (UserData*) user_data;
  int retval;

  /* NOTE: The order in which Advection and Reaction are
           called is critical here. Advection must be
           computed first. */
  retval = Advection(t, y, ydot, user_data);
  if (check_retval((void *)&retval, "Advection", 1, udata->myid)) return(-1);

  retval = Reaction(t, y, ydot, user_data);
  if (check_retval((void *)&retval, "Reaction", 1, udata->myid)) return(-1);

  /* return success */
  return(0);
}

/* Compute the residual F(t,y,y') = ydot - h(t,y) = 0. */
static int AdvectionReactionResidual(realtype t, N_Vector y, N_Vector ydot,
                                     N_Vector F, void *user_data)
{
  /* access problem data */
  UserData* udata = (UserData*) user_data;
  int retval;

  /* NOTE: The order in which Advection and Reaction are
           called is critical here. Advection must be
           computed first. */
  retval = Advection(t, y, F, user_data); /* F = -c y_x */
  if (check_retval((void *)&retval, "Advection", 1, udata->myid)) return(-1);

  retval = Reaction(t, y, F, user_data);  /* F = -c y_x + g(t,y) */
  if (check_retval((void *)&retval, "Reaction", 1, udata->myid)) return(-1);

  /* F = ydot - h(t,y) = ydot + c y_x - g(t,y) */
  N_VLinearSum(1.0, ydot, -1.0, F, F);

  /* return success */
  return(0);
}

/* --------------------------------------------------------------
 * Linear system and Jacobian functions
 * --------------------------------------------------------------*/

/* Solve the linear systems Ax = b where A = I - gamma*dg/dy.
   When using a fully implicit method, we are approximating
   dh/dy as dg/dy. */
static int SolveReactionLinSys(N_Vector y, N_Vector x, N_Vector b,
                               realtype gamma, raja_xyz_tuple blocks,
                               UserData* udata)
{
  /* shortcuts */
  int       dof, nxl, nyl, nzl;
  realtype  k2, k3, k4, k6;

  /* set shortcuts */
  dof = udata->grid->dof;
  nxl = udata->grid->nxl;
  nyl = udata->grid->nyl;
  nzl = udata->grid->nzl;
  k2  = udata->k2;
  k3  = udata->k3;
  k4  = udata->k4;
  k6  = udata->k6;
  
  /* create views of the data */
  RAJA::View<realtype, RAJA::Layout<NDIMS+1> > Yview(GetVecData(y),
                                                     nxl, nyl, nzl, dof);
  RAJA::View<realtype, RAJA::Layout<NDIMS+1> > Bview(GetVecData(b),
                                                     nxl, nyl, nzl, dof);
  RAJA::View<realtype, RAJA::Layout<NDIMS+1> > Xview(GetVecData(x),
                                                     nxl, nyl, nzl, dof);

  RAJA::kernel<XYZ_KERNEL_POL>(blocks,
    [=] DEVICE_FUNC (int i, int j, int k) {

    /* and the corresponding vectors */
    realtype *b = &(Bview(i,j,k,0));
    realtype *x = &(Xview(i,j,k,0));

    /* shortcuts to u, v, w for the block */
    realtype u = Yview(i,j,k,0);
    realtype v = Yview(i,j,k,1);
    realtype w = Yview(i,j,k,2);

    realtype A0, A1, A2, A3, A4, A5, A6, A7, A8;

    //
    // compute J = dg/dy
    //

    /* 1st row: u, v, w */
    A0 = -k2 * w + 2.0 * k3 * u * v - k4;
    A1 =  k3 * u * u;
    A2 = -k2 * u;

    /* 2nd row: u, v, w */
    A3 =  k2 * w - 2.0 * k3 * u * v;
    A4 = -k3 * u * u;
    A5 =  k2 * u;

    /* 3rd row: u, v, w */
    A6 = -k2 * w;
    A7 =  0.0;
    A8 = -k2 * u - k6;

    //
    // compute A = I - gamma*J
    //

    A0 = 1. - (gamma * A0);
    A1 = -gamma * A1;
    A2 = -gamma * A2;
    A3 = -gamma * A3;
    A4 = 1. - (gamma * A4);
    A5 = -gamma * A5;
    A6 = -gamma * A6;
    A7 = -gamma * A7;
    A8 = 1. - (gamma * A8);

    //
    // compute x = A^{-1}b
    //

    realtype scratch_0 = A4*A8;
    realtype scratch_1 = A1*A5;
    realtype scratch_2 = A2*A7;
    realtype scratch_3 = A5*A7;
    realtype scratch_4 = A1*A8;
    realtype scratch_5 = A2*A4;
    realtype scratch_6 = 1.0/(A0*scratch_0 - A0*scratch_3 + A3*scratch_2 - A3*scratch_4 + A6*scratch_1 - A6*scratch_5);
    realtype scratch_7 = A2*A3;
    realtype scratch_8 = A6*b[0];
    realtype scratch_9 = A2*A6;
    realtype scratch_10 = A3*b[0];
    realtype scratch_11 = 1.0/A0;
    realtype scratch_12 = A1*scratch_11;
    realtype scratch_13 = (-A6*scratch_12 + A7)/(-A3*scratch_12 + A4);

    x[0] = scratch_6*(b[0]*scratch_0 - b[0]*scratch_3 + b[1]*scratch_2 - b[1]*scratch_4 + b[2]*scratch_1 - b[2]*scratch_5);
    x[1] = scratch_6*(-A0*A5*b[2] + A0*A8*b[1] + A5*scratch_8 - A8*scratch_10 - b[1]*scratch_9 + b[2]*scratch_7);
    x[2] = (-b[2] + scratch_11*scratch_8 + scratch_13*(b[1] - scratch_10*scratch_11))/(-A8 + scratch_11*scratch_9 + scratch_13*(A5 - scratch_11*scratch_7));
  });

  return(0);
}

/* Solve the linear systems Ax = b where A = -dg/dy + gamma.
   We are approximating dh/dy as dg/dy. */
static int SolveReactionLinSysRes(N_Vector y, N_Vector x, N_Vector b,
                                  realtype gamma, raja_xyz_tuple blocks,
                                  UserData* udata)
{
  /* shortcuts */
  int       dof, nxl, nyl, nzl;
  realtype  k2, k3, k4, k6;

  /* set shortcuts */
  dof = udata->grid->dof;
  nxl = udata->grid->nxl;
  nyl = udata->grid->nyl;
  nzl = udata->grid->nzl;
  k2    = udata->k2;
  k3    = udata->k3;
  k4    = udata->k4;
  k6    = udata->k6;

  /* create views of the data */
  RAJA::View<realtype, RAJA::Layout<NDIMS+1> > Yview(GetVecData(y),
                                                     nxl, nyl, nzl, dof);
  RAJA::View<realtype, RAJA::Layout<NDIMS+1> > Bview(GetVecData(b),
                                                     nxl, nyl, nzl, dof);
  RAJA::View<realtype, RAJA::Layout<NDIMS+1> > Xview(GetVecData(x),
                                                     nxl, nyl, nzl, dof);

  RAJA::kernel<XYZ_KERNEL_POL>(blocks,
    [=] DEVICE_FUNC (int i, int j, int k) {

    /* and the corresponding vectors */
    realtype *b = &(Bview(i,j,k,0));
    realtype *x = &(Xview(i,j,k,0));

    /* shortcuts to u, v, w for the block */
    realtype u = Yview(i,j,k,0);
    realtype v = Yview(i,j,k,1);
    realtype w = Yview(i,j,k,2);

    realtype A0, A1, A2, A3, A4, A5, A6, A7, A8;

    //
    // compute dg/dy
    //

    /* 1st row: u, v, w */
    A0 = -k2 * w + 2.0 * k3 * u * v - k4;
    A1 =  k3 * u * u;
    A2 = -k2 * u;

    /* 2nd row: u, v, w */
    A3 =  k2 * w - 2.0 * k3 * u * v;
    A4 = -k3 * u * u;
    A5 =  k2 * u;

    /* 3rd row: u, v, w */
    A6 = -k2 * w;
    A7 =  0.0;
    A8 = -k2 * u - k6;

    //
    // compute A = -dg/dy + gamma*diag(df/dydot)
    // where diag(df/dydot) is approximated as
    // diag([udot, vdot, wdot])
    //

    A0 = -A0 + gamma;
    A1 = -A1;
    A2 = -A2;
    A3 = -A3;
    A4 = -A4 + gamma;
    A5 = -A5;
    A6 = -A6;
    A7 = -A7;
    A8 = -A8 + gamma;

    //
    // compute x = A^{-1}b
    //

    realtype scratch_0 = A4*A8;
    realtype scratch_1 = A1*A5;
    realtype scratch_2 = A2*A7;
    realtype scratch_3 = A5*A7;
    realtype scratch_4 = A1*A8;
    realtype scratch_5 = A2*A4;
    realtype scratch_6 = 1.0/(A0*scratch_0 - A0*scratch_3 + A3*scratch_2 - A3*scratch_4 + A6*scratch_1 - A6*scratch_5);
    realtype scratch_7 = A2*A3;
    realtype scratch_8 = A6*b[0];
    realtype scratch_9 = A2*A6;
    realtype scratch_10 = A3*b[0];
    realtype scratch_11 = 1.0/A0;
    realtype scratch_12 = A1*scratch_11;
    realtype scratch_13 = (-A6*scratch_12 + A7)/(-A3*scratch_12 + A4);

    x[0] = scratch_6*(b[0]*scratch_0 - b[0]*scratch_3 + b[1]*scratch_2 - b[1]*scratch_4 + b[2]*scratch_1 - b[2]*scratch_5);
    x[1] = scratch_6*(-A0*A5*b[2] + A0*A8*b[1] + A5*scratch_8 - A8*scratch_10 - b[1]*scratch_9 + b[2]*scratch_7);
    x[2] = (-b[2] + scratch_11*scratch_8 + scratch_13*(b[1] - scratch_10*scratch_11))/(-A8 + scratch_11*scratch_9 + scratch_13*(A5 - scratch_11*scratch_7));
  });

  return(0);
}


/* --------------------------------------------------------------
 * Preconditioner functions
 * --------------------------------------------------------------*/

/* Solves Pz = r where P = I - gamma * dg/dy */
static int PSolve(realtype t, N_Vector y, N_Vector ydot, N_Vector r,
                  N_Vector z, realtype gamma, realtype delta, int lr,

                  void *user_data)
{
  /* local variables */
  UserData* udata = (UserData*) user_data;
  int       retval;

  SUNDIALS_CXX_MARK_FUNCTION(udata->prof);

  /* solve the task-local linear system Pz = r */
  auto range = RAJA::make_tuple(RAJA::RangeSegment(0, udata->grid->nxl),
                                RAJA::RangeSegment(0, udata->grid->nyl),
                                RAJA::RangeSegment(0, udata->grid->nzl));
  retval = SolveReactionLinSys(y, z, r, gamma, range, udata);

  return(retval);
}

/* Solves Pz = r where P = -dg/dy + gamma */
static int PSolveRes(realtype t, N_Vector y, N_Vector ydot, N_Vector F,
                     N_Vector r, N_Vector z, realtype cj, realtype delta,
                     void *user_data)
{
  /* local variables */
  UserData* udata = (UserData*) user_data;
  int       retval;

  SUNDIALS_CXX_MARK_FUNCTION(udata->prof);

  /* solve the task-local linear system Pz = r */
  auto range = RAJA::make_tuple(RAJA::RangeSegment(0, udata->grid->nxl),
                                RAJA::RangeSegment(0, udata->grid->nyl),
                                RAJA::RangeSegment(0, udata->grid->nzl));
  retval = SolveReactionLinSysRes(y, z, r, cj, range, udata);

  return(retval);
}


#endif
