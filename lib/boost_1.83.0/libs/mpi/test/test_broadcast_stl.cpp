// Copyright (C) 2005, 2006 Douglas Gregor.

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

// A test of the broadcast() collective.
#include <algorithm>
#include <vector>
#include <map>

#include <boost/mpi/collectives/broadcast.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>

#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>

#include "mpi_test_utils.hpp"

namespace mpi = boost::mpi;

typedef std::vector<std::map<int, double> > sparse;

template<typename T>
int
broadcast_test(const mpi::communicator& comm, const T& bc_value,
               std::string const& kind, int root) {
  using boost::mpi::broadcast;
  int failed = 0;
  T value;
  if (comm.rank() == root) {
    value = bc_value;
    std::cout << "Broadcasting " << kind << " from root " << root << "...";
    std::cout.flush();
  }
  
  broadcast(comm, value, root);
  BOOST_MPI_CHECK(value == bc_value, failed);
  if (comm.rank() == root) {
    if (value == bc_value) {
      std::cout << "OK." << std::endl;
    } else {
      std::cout << "FAIL." << std::endl;
    }
  }
  comm.barrier();
  return failed;
}

template<typename T>
int
broadcast_test(const mpi::communicator& comm, const T& bc_value,
               std::string const& kind)
{
  int failed = 0;
  for (int root = 0; root < comm.size(); ++root) {
    BOOST_MPI_COUNT_FAILED(broadcast_test(comm, bc_value, kind, root), failed);
  }
  return failed;
}

int main()
{
  boost::mpi::environment env;
  int failed = 0;
  mpi::communicator comm;
  BOOST_MPI_CHECK(comm.size() > 1, failed);
  if (failed == 0) {
    sparse s;
    s.resize(2);
    s[0][12] = 0.12;
    s[1][13] = 1.13;
    BOOST_MPI_COUNT_FAILED(broadcast_test(comm, s, "sparse"), failed);
  }
  return failed;
}
