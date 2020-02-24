#ifdef STAN_LANG_MPI

#include <gtest/gtest.h>
#include <stan/math/prim/mat/fun/welford_var_estimator.hpp>
#include <stan/math/mpi/mpi_var_estimator.hpp>
#include <boost/mpi.hpp>
#include <fstream>

using Eigen::MatrixXd;
using Eigen::Matrix;
using std::vector;
using stan::math::mpi::Communicator;
using stan::math::mpi::Session;

TEST(mpi_variance_test, var_estimator) {
  const Communicator comm(MPI_COMM_STAN);

  const int num_draws = 3;
  const int num_params = 5;
  const int max_num_chains = 4;
  const int num_chains = comm.size();
  Eigen::VectorXd draws_pool(num_draws * num_params * max_num_chains);
  draws_pool <<
    -266.132,  -265.76 ,  -265.757,  -265.793,  -265.791,
    -261.793,  -265.598,  -265.354,  -267.131,  -265.039,
    -265.039,  -265.031,  -261.039,  -262.128,  -260.869,
    -265.869,  -266.309,  -265.897,  -265.727,  -265.958,
    -267.231,  -266.862,  -267.255,  -267.545,  -267.009,
    -265.8  ,  -265.551,  -266.254,  -265.394,  -264.825,
    -264.825,  -265.205,  -262.245,  -261.312,  -268.365,
    -264.253,  -264.514,  -264.413,  -264.413,  -264.413,
    -268.413,  -264.413,  -264.589,  -265.277,  -265.378,
    -265.69 ,  -265.69 ,  -264.972,  -264.972,  -269.972,
    -264.972,  -265.283,  -265.237,  -264.671,  -264.88 ,
    -265.099,  -266.919,  -265.878,  -269.653,  -264.653;

  Eigen::VectorXd all_draws(num_draws * num_params * num_chains);
  for (int i = 0; i < all_draws.size(); ++i) {
    all_draws(i) = draws_pool(i);
  }

  std::vector<Eigen::VectorXd> chain_draws(num_draws);
  for (int i = 0; i < num_draws; ++i) {
    int i_begin = comm.rank() * num_draws * num_params + i * num_params;
    chain_draws[i] = Eigen::VectorXd::Map(&all_draws(i_begin), num_params);
  }

  stan::math::mpi::mpi_var_estimator estimator(num_params);
  for (int i = 0; i < num_draws; ++i) {
    estimator.add_sample(chain_draws[i]);
  }

  EXPECT_EQ(num_draws * comm.size(), estimator.num_samples(comm)[1]);

  stan::math::welford_var_estimator welford(num_params);
  for (int i = 0; i < num_draws * num_chains; ++i) {
    Eigen::VectorXd var(Eigen::VectorXd::Map(&all_draws(i * num_params), num_params));
    welford.add_sample(var);
  }

  Eigen::VectorXd m, m0;
  Eigen::VectorXd v, v0;
  estimator.sample_mean(m, comm);
  welford.sample_mean(m0);
  for (int i = 0; i < num_params; ++i) {
    EXPECT_FLOAT_EQ(m(i), m0(i));
  }

  estimator.sample_variance(v, comm);
  welford.sample_variance(v0);
  for (int i = 0; i < num_params; ++i) {
    EXPECT_FLOAT_EQ(v(i), v0(i));
  }
}

#endif
