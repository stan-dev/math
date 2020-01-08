#ifdef STAN_LANG_MPI

#include <gtest/gtest.h>
#include <stan/math/mpi/envionment.hpp>
#include <boost/mpi.hpp>
#include <fstream>

using Eigen::MatrixXd;
using Eigen::Matrix;
using std::vector;
using stan::math::mpi::Communicator;
using stan::math::mpi::Session;

TEST(mpi_comm_test, mpi_inter_intra_comms) {
  const Communicator world_comm(MPI_COMM_STAN);
  const Communicator& inter_comm(Session::inter_chain_comm(3));
  const Communicator& intra_comm(Session::intra_chain_comm(3));
  if (world_comm.size() == 3) {
    switch (world_comm.rank()) {
    case 0:
      EXPECT_EQ(inter_comm.rank(), 0);
      EXPECT_EQ(intra_comm.rank(), 0);
      break;
    case 1:
      EXPECT_EQ(inter_comm.rank(), 1);
      EXPECT_EQ(intra_comm.rank(), 0);
      break;
    case 2:
      EXPECT_EQ(inter_comm.rank(), 2);
      EXPECT_EQ(intra_comm.rank(), 0);
      break;
    }
  } else if (world_comm.size() == 4) {
    switch (world_comm.rank()) {
    case 0:
      EXPECT_EQ(inter_comm.rank(), 0);
      EXPECT_EQ(intra_comm.rank(), 0);
      break;
    case 1:
      EXPECT_EQ(inter_comm.rank(), -1);
      EXPECT_EQ(intra_comm.rank(), 1);
      break;
    case 2:
      EXPECT_EQ(inter_comm.rank(), 1);
      EXPECT_EQ(intra_comm.rank(), 0);
      break;
    case 3:
      EXPECT_EQ(inter_comm.rank(), 2);
      EXPECT_EQ(intra_comm.rank(), 0);
      break;
    }
  } else if (world_comm.size() == 5) {
    switch (world_comm.rank()) {
    case 0:
      EXPECT_EQ(inter_comm.rank(), 0);
      EXPECT_EQ(intra_comm.rank(), 0);
      break;
    case 1:
      EXPECT_EQ(inter_comm.rank(), -1);
      EXPECT_EQ(intra_comm.rank(), 1);
      break;
    case 2:
      EXPECT_EQ(inter_comm.rank(), 1);
      EXPECT_EQ(intra_comm.rank(), 0);
      break;
    case 3:
      EXPECT_EQ(inter_comm.rank(), -1);
      EXPECT_EQ(intra_comm.rank(), 1);
      break;
    case 4:
      EXPECT_EQ(inter_comm.rank(), 2);
      EXPECT_EQ(intra_comm.rank(), 0);
      break;
    }
  } else if (world_comm.size() == 6) {
    switch (world_comm.rank()) {
    case 0:
      EXPECT_EQ(inter_comm.rank(), 0);
      EXPECT_EQ(intra_comm.rank(), 0);
      break;
    case 1:
      EXPECT_EQ(inter_comm.rank(), -1);
      EXPECT_EQ(intra_comm.rank(), 1);
      break;
    case 2:
      EXPECT_EQ(inter_comm.rank(), 1);
      EXPECT_EQ(intra_comm.rank(), 0);
      break;
    case 3:
      EXPECT_EQ(inter_comm.rank(), -1);
      EXPECT_EQ(intra_comm.rank(), 1);
      break;
    case 4:
      EXPECT_EQ(inter_comm.rank(), 2);
      EXPECT_EQ(intra_comm.rank(), 0);
      break;
    case 5:
      EXPECT_EQ(inter_comm.rank(), -1);
      EXPECT_EQ(intra_comm.rank(), 1);
      break;
    }
  } else if (world_comm.size() == 7) {
    switch (world_comm.rank()) {
    case 0:
      EXPECT_EQ(inter_comm.rank(), 0);
      EXPECT_EQ(intra_comm.rank(), 0);
      break;
    case 1:
      EXPECT_EQ(inter_comm.rank(), -1);
      EXPECT_EQ(intra_comm.rank(), 1);
      break;
    case 2:
      EXPECT_EQ(inter_comm.rank(), -1);
      EXPECT_EQ(intra_comm.rank(), 2);
      break;
    case 3:
      EXPECT_EQ(inter_comm.rank(), 1);
      EXPECT_EQ(intra_comm.rank(), 0);
      break;
    case 4:
      EXPECT_EQ(inter_comm.rank(), -1);
      EXPECT_EQ(intra_comm.rank(), 1);
      break;
    case 5:
      EXPECT_EQ(inter_comm.rank(), 2);
      EXPECT_EQ(intra_comm.rank(), 0);
      break;
    case 6:
      EXPECT_EQ(inter_comm.rank(), -1);
      EXPECT_EQ(intra_comm.rank(), 1);
      break;
    }
  }
}

TEST(mpi_comm_test, is_in_inter_comm) {
  const Communicator world_comm(MPI_COMM_STAN);
  const int n_chains = 3;
  if (world_comm.size() == 3) {
    switch (world_comm.rank()) {
    case 0:
      EXPECT_TRUE(Session::is_in_inter_chain_comm(n_chains));
      break;
    case 1:
      EXPECT_TRUE(Session::is_in_inter_chain_comm(n_chains));
      break;
    case 2:
      EXPECT_TRUE(Session::is_in_inter_chain_comm(n_chains));
      break;
    }
  } else if (world_comm.size() == 4) {
    switch (world_comm.rank()) {
    case 0:
      EXPECT_TRUE(Session::is_in_inter_chain_comm(n_chains));
      break;
    case 1:
      EXPECT_FALSE(Session::is_in_inter_chain_comm(n_chains));
      break;
    case 2:
      EXPECT_TRUE(Session::is_in_inter_chain_comm(n_chains));
      break;
    case 3:
      EXPECT_TRUE(Session::is_in_inter_chain_comm(n_chains));
      break;
    }
  } else if (world_comm.size() == 5) {
    switch (world_comm.rank()) {
    case 0:
      EXPECT_TRUE(Session::is_in_inter_chain_comm(n_chains));
      break;
    case 1:
      EXPECT_FALSE(Session::is_in_inter_chain_comm(n_chains));
      break;
    case 2:
      EXPECT_TRUE(Session::is_in_inter_chain_comm(n_chains));
      break;
    case 3:
      EXPECT_FALSE(Session::is_in_inter_chain_comm(n_chains));
      break;
    case 4:
      EXPECT_TRUE(Session::is_in_inter_chain_comm(n_chains));
      break;
    }
  } else if (world_comm.size() == 6) {
    switch (world_comm.rank()) {
    case 0:
      EXPECT_TRUE(Session::is_in_inter_chain_comm(n_chains));
      break;
    case 1:
      EXPECT_FALSE(Session::is_in_inter_chain_comm(n_chains));
      break;
    case 2:
      EXPECT_TRUE(Session::is_in_inter_chain_comm(n_chains));
      break;
    case 3:
      EXPECT_FALSE(Session::is_in_inter_chain_comm(n_chains));
      break;
    case 4:
      EXPECT_TRUE(Session::is_in_inter_chain_comm(n_chains));
      break;
    case 5:
      EXPECT_FALSE(Session::is_in_inter_chain_comm(n_chains));
      break;
    }
  } else if (world_comm.size() == 7) {
    switch (world_comm.rank()) {
    case 0:
      EXPECT_TRUE(Session::is_in_inter_chain_comm(n_chains));
      break;
    case 1:
      EXPECT_FALSE(Session::is_in_inter_chain_comm(n_chains));
      break;
    case 2:
      EXPECT_FALSE(Session::is_in_inter_chain_comm(n_chains));
      break;
    case 3:
      EXPECT_TRUE(Session::is_in_inter_chain_comm(n_chains));
      break;
    case 4:
      EXPECT_FALSE(Session::is_in_inter_chain_comm(n_chains));
      break;
    case 5:
      EXPECT_TRUE(Session::is_in_inter_chain_comm(n_chains));
      break;
    case 6:
      EXPECT_FALSE(Session::is_in_inter_chain_comm(n_chains));
      break;
    }
  }
}

#endif
