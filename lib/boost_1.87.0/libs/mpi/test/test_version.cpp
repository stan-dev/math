// Copyright (C) 2013 Alain Miniussi <alain.miniussi@oca.eu>

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

// test mpi version

#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <iostream>

namespace mpi = boost::mpi;
// return number of failures
int
test_version(mpi::communicator const& comm) {
#if defined(MPI_VERSION)
  int mpi_version    = MPI_VERSION;
  int mpi_subversion = MPI_SUBVERSION;
#else
  int mpi_version = 0;
  int mpi_subversion = 0;
#endif
  
  std::pair<int,int> version = mpi::environment::version();
  if (comm.rank() == 0) {
    std::cout << "MPI Version: " << version.first << ',' << version.second << '\n';
  }
  int failed = 0;
  if (version.first != mpi_version) {
    ++failed;
  }
  if (version.second != mpi_subversion) {
    ++failed;
  }
  std::string lib_version = mpi::environment::library_version();
#if (3 <= MPI_VERSION)
  if (lib_version.size() == 0) {
    ++failed;
  }
#else
  if (lib_version.size() != 0) {
    ++failed;
  }
#endif
  return failed;
}

std::string
yesno(bool b) {
  return b ? std::string("yes") : std::string("no");
}

// not really a check
void
report_features(mpi::communicator const& comm) {
  if (comm.rank() == 0) {
    std::cout << "Assuming working MPI_Improbe:" <<
#if defined(BOOST_MPI_USE_IMPROBE)
      "yes" << '\n';
#else
      "no"  << '\n';
#endif
  }
}

int
main() {
  mpi::environment env;
  mpi::communicator world;

  int failed = 0;
  failed += test_version(world);
  report_features(world);
  return failed;
}
