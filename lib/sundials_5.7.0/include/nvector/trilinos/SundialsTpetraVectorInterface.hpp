/* -----------------------------------------------------------------
 * Programmer(s): Slaven Peles @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------*/

#ifndef _TPETRA_SUNDIALS_INTERFACE_HPP_
#define _TPETRA_SUNDIALS_INTERFACE_HPP_

#include <Tpetra_Vector.hpp>
#include <nvector/nvector_trilinos.h>

namespace Sundials
{

  struct TpetraVectorInterface : public _N_VectorContent_Trilinos
  {
    // Typedef of Tpetra vector class to be used with SUNDIALS
    typedef Tpetra::Vector<realtype, int, sunindextype> vector_type;

    TpetraVectorInterface(Teuchos::RCP<vector_type> rcpvec)
    {
      rcpvec_ = rcpvec;
    }

    ~TpetraVectorInterface() = default;

    Teuchos::RCP<vector_type> rcpvec_;
  };


} // namespace Sundials

inline Teuchos::RCP<Sundials::TpetraVectorInterface::vector_type> N_VGetVector_Trilinos(N_Vector v)
{
  Sundials::TpetraVectorInterface* iface =
  reinterpret_cast<Sundials::TpetraVectorInterface*>(v->content);

  return iface->rcpvec_;
}

/*
 * -----------------------------------------------------------------
 * Function : N_VMake_Trilinos
 * -----------------------------------------------------------------
 * This function attaches N_Vector functions to a Tpetra vector.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector
N_VMake_Trilinos(Teuchos::RCP<Sundials::TpetraVectorInterface::vector_type> v);



#endif // _TPETRA_SUNDIALS_INTERFACE_HPP_
