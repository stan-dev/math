#ifndef TEST_UNIT_MATH_PRIM_MAT_FUNCTOR_MPI_TEST_ENV_HPP
#define TEST_UNIT_MATH_PRIM_MAT_FUNCTOR_MPI_TEST_ENV_HPP

#include <gtest/gtest.h>

// sets up the MPI environment. All tests have to be skipped on
// non-root nodes with a
// if(rank != 0) return;
// line at the test start

// moreover, the MPI tests have to be setup with the MPI_TEST_F macro
// which will disable all MPI tests whenever no STAN_HAS_MPI is
// defined

#ifdef STAN_MPI

#define MPI_TEST_F(test_fixture, test_name) TEST_F(test_fixture, test_name)

#define MPI_TEST(test_fixture, test_name) TEST(test_fixture, test_name)

#include <stan/math/prim/arr/functor/mpi_cluster.hpp>

stan::math::mpi_cluster* cluster;
std::size_t rank;
std::size_t world_size;

class MPIEnvironment : public ::testing::Environment {
 public:
  virtual void SetUp() {
    cluster = new stan::math::mpi_cluster;
    boost::mpi::communicator world;
    rank = world.rank();
    world_size = world.size();
    ::testing::TestEventListeners& listeners
        = ::testing::UnitTest::GetInstance()->listeners();
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
::testing::Environment* const mpi_env
    = ::testing::AddGlobalTestEnvironment(new MPIEnvironment);

#else

#define MPI_TEST_F(test_fixture, test_name) \
  TEST_F(test_fixture, DISABLED_##test_name)

#define MPI_TEST(test_fixture, test_name) \
  TEST(test_fixture, DISABLED_##test_name)

std::size_t rank = 0;
std::size_t world_size = 0;

#endif

#endif
