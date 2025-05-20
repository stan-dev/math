// Copyright (C) 2005, 2006 Douglas Gregor.

// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

// A test of the communicator that passes data around a ring and
// verifies that the same data makes it all the way. Should test all
// of the various kinds of data that can be sent (primitive types, POD
// types, serializable objects, etc.)
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>
#include <algorithm>
#include "gps_position.hpp"
#include <boost/serialization/string.hpp>
#include <boost/serialization/list.hpp>
//#include "debugger.cpp"

#include "mpi_test_utils.hpp"

using boost::mpi::communicator;
using boost::mpi::status;

template<typename T>
int
ring_test(const communicator& comm, const T& pass_value, const char* kind,
          int root = 0)
{
  int failed  = 0;
  T transferred_value;

  int rank = comm.rank();
  int size = comm.size();

  if (rank == root) {
    std::cout << "Passing " << kind << " around a ring from root " << root
              << "...";
    comm.send((rank + 1) % size, 0, pass_value);
    comm.recv((rank + size - 1) % size, 0, transferred_value);
    BOOST_MPI_CHECK(transferred_value == pass_value, failed);
    if (transferred_value == pass_value) std::cout << " OK." << std::endl;
  } else {
    comm.recv((rank + size - 1) % size, 0, transferred_value);
    BOOST_MPI_CHECK(transferred_value == pass_value, failed);
    comm.send((rank + 1) % size, 0, transferred_value);
  }
  
  (comm.barrier)();
  return failed;
}


template<typename T>
int
ring_array_test(const communicator& comm, const T* pass_values,
                int n, const char* kind, int root = 0)
{
  int failed = 0;
  T* transferred_values = new T[n];
  int rank = comm.rank();
  int size = comm.size();

  if (rank == root) {

    std::cout << "Passing " << kind << " array around a ring from root "
              << root  << "...";
    comm.send((rank + 1) % size, 0, pass_values, n);
    comm.recv((rank + size - 1) % size, 0, transferred_values, n);
    bool okay = std::equal(pass_values, pass_values + n,
                           transferred_values);
    BOOST_MPI_CHECK(okay, failed);
    if (okay) std::cout << " OK." << std::endl;
  } else {
    status stat = comm.probe(boost::mpi::any_source, 0);
    boost::optional<int> num_values = stat.template count<T>();
    if (boost::mpi::is_mpi_datatype<T>()) {
      BOOST_MPI_CHECK(num_values && *num_values == n, failed);
    } else {
      BOOST_MPI_CHECK(!num_values || *num_values == n, failed);     
    }
    comm.recv(stat.source(), 0, transferred_values, n);
    BOOST_MPI_CHECK(std::equal(pass_values, pass_values + n,
                         transferred_values), failed);
    comm.send((rank + 1) % size, 0, transferred_values, n);
  }
  (comm.barrier)();
  delete [] transferred_values;
  return failed;
}

enum color_t {red, green, blue};
BOOST_IS_MPI_DATATYPE(color_t)

int main()
{
  boost::mpi::environment env;
  communicator comm;
  int failed = 0;
  BOOST_MPI_CHECK(comm.size() > 1, failed);
  if (failed == 0) {
    // Check transfer of individual objects
    BOOST_MPI_COUNT_FAILED(ring_test(comm, 17, "integers", 0), failed);
    BOOST_MPI_COUNT_FAILED(ring_test(comm, 17, "integers", 1), failed);
    BOOST_MPI_COUNT_FAILED(ring_test(comm, red, "enums", 1), failed);
    BOOST_MPI_COUNT_FAILED(ring_test(comm, red, "enums", 1), failed);
    BOOST_MPI_COUNT_FAILED(ring_test(comm, gps_position(39,16,20.2799), "GPS positions", 0), failed);
    BOOST_MPI_COUNT_FAILED(ring_test(comm, gps_position(26,25,30.0), "GPS positions", 1), failed);
    BOOST_MPI_COUNT_FAILED(ring_test(comm, std::string("Rosie"), "string", 0), failed);
    
    std::list<std::string> strings;
    strings.push_back("Hello");
    strings.push_back("MPI");
    strings.push_back("World");
    BOOST_MPI_COUNT_FAILED(ring_test(comm, strings, "list of strings", 1), failed);
    
    // Check transfer of arrays
    int int_array[2] = { 17, 42 };
    BOOST_MPI_COUNT_FAILED(ring_array_test(comm, int_array, 2, "integer", 1), failed);
    gps_position gps_position_array[2] = {
      gps_position(39,16,20.2799),
      gps_position(26,25,30.0)
    };
    BOOST_MPI_COUNT_FAILED(ring_array_test(comm, gps_position_array, 2, "GPS position", 1), failed);
    
    std::string string_array[3] = { "Hello", "MPI", "World" };
    BOOST_MPI_COUNT_FAILED(ring_array_test(comm, string_array, 3, "string", 0), failed);
    BOOST_MPI_COUNT_FAILED(ring_array_test(comm, string_array, 3, "string", 1), failed);
  }
  return failed;
}
