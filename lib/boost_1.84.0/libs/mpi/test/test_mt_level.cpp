// Copyright (C) 2013 Alain Miniussi <alain.miniussi@oca.eu>

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

// test threading::level operations

#include <boost/mpi/environment.hpp>
#include <iostream>
#include <sstream>

#include "mpi_test_utils.hpp"

namespace mpi = boost::mpi;

int 
test_threading_level_io(mpi::threading::level orig) {
  int failed = 0;
  std::ostringstream out;
  namespace mt = boost::mpi::threading;
  mt::level printed = mt::level(-1);

  out << orig;
  BOOST_MPI_CHECK(out.good(), failed);
  std::string orig_str(out.str());
  std::cout << "orig string:" << orig_str << '\n';
  std::istringstream in(orig_str);
  in >> printed;
  BOOST_MPI_CHECK(!in.bad(), failed);
  std::cout << "orig: " << orig << ", printed: " << printed << std::endl;
  BOOST_MPI_CHECK(orig == printed, failed);
  return failed;
}

int
test_threading_levels_io() {
  namespace mt = boost::mpi::threading;
  int failed = 0;
  BOOST_MPI_COUNT_FAILED(test_threading_level_io(mt::single), failed);
  BOOST_MPI_COUNT_FAILED(test_threading_level_io(mt::funneled), failed);
  BOOST_MPI_COUNT_FAILED(test_threading_level_io(mt::serialized), failed);
  BOOST_MPI_COUNT_FAILED(test_threading_level_io(mt::multiple), failed);
  return failed;
}

int
test_threading_level_cmp() {
  namespace mt = boost::mpi::threading;
  int failed = 0;
  BOOST_MPI_CHECK(mt::single == mt::single, failed);
  BOOST_MPI_CHECK(mt::funneled == mt::funneled, failed);
  BOOST_MPI_CHECK(mt::serialized == mt::serialized, failed);
  BOOST_MPI_CHECK(mt::multiple == mt::multiple, failed);
  
  BOOST_MPI_CHECK(mt::single != mt::funneled, failed);
  BOOST_MPI_CHECK(mt::single != mt::serialized, failed);
  BOOST_MPI_CHECK(mt::single != mt::multiple, failed);
  
  BOOST_MPI_CHECK(mt::funneled != mt::single, failed);
  BOOST_MPI_CHECK(mt::funneled != mt::serialized, failed);
  BOOST_MPI_CHECK(mt::funneled != mt::multiple, failed);
  
  BOOST_MPI_CHECK(mt::serialized != mt::single, failed);
  BOOST_MPI_CHECK(mt::serialized != mt::funneled, failed);
  BOOST_MPI_CHECK(mt::serialized != mt::multiple, failed);
  
  BOOST_MPI_CHECK(mt::multiple != mt::single, failed);
  BOOST_MPI_CHECK(mt::multiple != mt::funneled, failed);
  BOOST_MPI_CHECK(mt::multiple != mt::serialized, failed);
  
  BOOST_MPI_CHECK(mt::single < mt::funneled, failed);
  BOOST_MPI_CHECK(mt::funneled > mt::single, failed);
  BOOST_MPI_CHECK(mt::single < mt::serialized, failed);
  BOOST_MPI_CHECK(mt::serialized > mt::single, failed);
  BOOST_MPI_CHECK(mt::single < mt::multiple, failed);
  BOOST_MPI_CHECK(mt::multiple > mt::single, failed);
  
  BOOST_MPI_CHECK(mt::funneled < mt::serialized, failed);
  BOOST_MPI_CHECK(mt::serialized > mt::funneled, failed);
  BOOST_MPI_CHECK(mt::funneled < mt::multiple, failed);
  BOOST_MPI_CHECK(mt::multiple > mt::funneled, failed);
  
  BOOST_MPI_CHECK(mt::serialized < mt::multiple, failed);
  BOOST_MPI_CHECK(mt::multiple > mt::serialized, failed);
  
  BOOST_MPI_CHECK(mt::single <= mt::single, failed);
  BOOST_MPI_CHECK(mt::single <= mt::funneled, failed);
  BOOST_MPI_CHECK(mt::funneled >= mt::single, failed);
  BOOST_MPI_CHECK(mt::single <= mt::serialized, failed);
  BOOST_MPI_CHECK(mt::serialized >= mt::single, failed);
  BOOST_MPI_CHECK(mt::single <= mt::multiple, failed);
  BOOST_MPI_CHECK(mt::multiple >= mt::single, failed);
  
  BOOST_MPI_CHECK(mt::funneled <= mt::funneled, failed);
  BOOST_MPI_CHECK(mt::funneled <= mt::serialized, failed);
  BOOST_MPI_CHECK(mt::serialized >= mt::funneled, failed);
  BOOST_MPI_CHECK(mt::funneled <= mt::multiple, failed);
  BOOST_MPI_CHECK(mt::multiple >= mt::funneled, failed);
  
  BOOST_MPI_CHECK(mt::serialized <= mt::serialized, failed);
  BOOST_MPI_CHECK(mt::serialized <= mt::multiple, failed);
  BOOST_MPI_CHECK(mt::multiple >= mt::serialized, failed);
  
  BOOST_MPI_CHECK(mt::multiple <= mt::multiple, failed);
  return failed;
}
    
int main()
{
  int failed = 0;
  BOOST_MPI_COUNT_FAILED(test_threading_levels_io(), failed);
  BOOST_MPI_COUNT_FAILED(test_threading_level_cmp(), failed);
  return failed;
}
