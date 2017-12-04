#pragma once

// sets up the MPI environment. All tests have to be skipped on
// non-root nodes with a
// if(rank != 0) return
// line at the test start

#include <gtest/gtest.h>

#include <stan/math/prim/arr/functor/mpi_cluster.hpp>

stan::math::mpi_cluster* cluster;
std::size_t rank;

class MPIEnvironment : public ::testing::Environment
{
public:
  
  virtual void SetUp() {
    cluster = new stan::math::mpi_cluster;
    boost::mpi::communicator world;
    rank = world.rank();
    ::testing::TestEventListeners& listeners =
        ::testing::UnitTest::GetInstance()->listeners();
    if (rank != 0) {
      delete listeners.Release(listeners.default_result_printer());
    }
    cluster->listen();
  }
  virtual void TearDown() {
    cluster->stop_listen();
    delete cluster;
  }

  virtual ~MPIEnvironment() {}
};


// register MPI global
 ::testing::Environment* const mpi_env = ::testing::AddGlobalTestEnvironment(new MPIEnvironment);

