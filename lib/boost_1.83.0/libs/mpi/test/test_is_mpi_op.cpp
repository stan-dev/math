// Copyright (C) 2005-2006 Douglas Gregor <doug.gregor@gmail.com>

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

// A test of the is_mpi_op functionality.
#include <boost/mpi/operations.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/type_traits/is_base_and_derived.hpp>

#include "mpi_test_utils.hpp"

using namespace boost::mpi;
using namespace std;
using boost::is_base_and_derived;

int main()
{
  boost::mpi::environment env;
  int failed = 0;
  // Check each predefined MPI_Op type that we support directly.
  BOOST_MPI_CHECK((is_mpi_op<minimum<float>, float>::op() == MPI_MIN), failed);
  BOOST_MPI_CHECK((is_mpi_op<plus<double>, double>::op() == MPI_SUM), failed);
  BOOST_MPI_CHECK((is_mpi_op<multiplies<long>, long>::op() == MPI_PROD), failed);
  BOOST_MPI_CHECK((is_mpi_op<logical_and<int>, int>::op() == MPI_LAND), failed);
  BOOST_MPI_CHECK((is_mpi_op<bitwise_and<int>, int>::op() == MPI_BAND), failed);
  BOOST_MPI_CHECK((is_mpi_op<logical_or<int>, int>::op() == MPI_LOR), failed);
  BOOST_MPI_CHECK((is_mpi_op<bitwise_or<int>, int>::op() == MPI_BOR), failed);
  BOOST_MPI_CHECK((is_mpi_op<logical_xor<int>, int>::op() == MPI_LXOR), failed);
  BOOST_MPI_CHECK((is_mpi_op<bitwise_xor<int>, int>::op() == MPI_BXOR), failed);
  return failed;
}
