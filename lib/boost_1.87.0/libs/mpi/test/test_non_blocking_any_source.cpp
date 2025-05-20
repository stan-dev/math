// Copyright (C) 2018 Steffen Hirschmann

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

// Test any_source on serialized data
#include <vector>
#include <iostream>
#include <iterator>
#include <boost/mpi.hpp>
#include <boost/serialization/vector.hpp>

#include "mpi_test_utils.hpp"

namespace mpi = boost::mpi;

std::string ok(bool b) {
  return b ? "ok" : "ko";
}

int main()
{
  mpi::environment env;
  mpi::communicator world;
  int failed = 0;
  int rank = world.rank();
#if BOOST_MPI_VERSION < 3
  if (rank == 0) {
    std::cout << "\nExpected failure with MPI standard < 3 ("
              << BOOST_MPI_VERSION << "." << BOOST_MPI_SUBVERSION
              << " detected)\n\n";
  }
  return;
#endif
  if (rank == 0) {
    std::vector<boost::mpi::request> req;
    std::vector<std::vector<int> > data(world.size() - 1);
    for (int i = 1; i < world.size(); ++i) {
      req.push_back(world.irecv(mpi::any_source, 0, data[i - 1]));
    }
    boost::mpi::wait_all(req.begin(), req.end());
    std::vector<bool> check(world.size()-1, false);
    for (int i = 0; i < world.size() - 1; ++i) {
      std::cout << "Process 0 received:" << std::endl;
      std::copy(data[i].begin(), data[i].end(), std::ostream_iterator<int>(std::cout, " "));
      std::cout << std::endl;
      int idx = data[i].size();
      BOOST_MPI_CHECK(std::equal_range(data[i].begin(), data[i].end(), idx)
                == std::make_pair(data[i].begin(), data[i].end()), failed);
      check[idx-1] = true;
    }
    for(int i = 0; i < world.size() - 1; ++i) {
      std::cout << "Received from " << i+1 << " is " << ok(check[i]) << '\n';
    }
    BOOST_MPI_CHECK(std::equal_range(check.begin(), check.end(), true)
              == std::make_pair(check.begin(), check.end()), failed); 
  } else {
    std::vector<int> vec(rank, rank);
    mpi::request req = world.isend(0, 0, vec);
    req.wait();
  }
  return failed;
}
