// Copyright (C) 2017 Alain Miniussi & Vincent Chabannes

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include <iostream>
#include <sstream>
#include <vector>

#include <boost/mpi.hpp>
#include <boost/mpi/nonblocking.hpp>
#include <boost/serialization/string.hpp>

#include "mpi_test_utils.hpp"

namespace mpi = boost::mpi;

int main()
{
  mpi::environment  env;
  mpi::communicator comm;
  int failed = 0;
  int rank = comm.rank();
  int const sz = 10;
  std::vector<int> data;
  std::vector< mpi::request> reqs;
  if ( rank == 0 ) {
    for ( int i=0; i<sz; ++i ) {
      data.push_back( i );
    }
    reqs.push_back( comm.isend(1, 0, data) );
  } else if ( rank == 1 ) {
    reqs.push_back( comm.irecv(0, 0, data) );
  }
  mpi::wait_all( reqs.begin(), reqs.end() );
  
  if ( rank == 1 ) {
    BOOST_MPI_CHECK(data.size() == sz, failed);
    for ( int i=0; i<sz; ++i ) {
      BOOST_MPI_CHECK(data[i] == i, failed);
    }
  }  
  return failed;
}
