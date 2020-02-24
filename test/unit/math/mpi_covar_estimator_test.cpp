#ifdef STAN_LANG_MPI

#include <gtest/gtest.h>
#include <stan/math/prim/mat/fun/welford_covar_estimator.hpp>
#include <stan/math/mpi/mpi_covar_estimator.hpp>
#include <boost/mpi.hpp>
#include <fstream>

using Eigen::MatrixXd;
using Eigen::Matrix;
using std::vector;
using stan::math::mpi::Communicator;
using stan::math::mpi::Session;

TEST(mpi_covariance_test, covar_using_all_draws) {
  const Communicator comm(MPI_COMM_STAN);

  const int num_draws = 3;
  const int num_params = 5;
  const int max_num_chains = 4;
  const int num_chains = comm.size();
  Eigen::VectorXd draws_pool(num_draws * num_params * max_num_chains);
  draws_pool <<
    -266.132,  -65.76 ,  265.757,  -2.793,  -265.791,
    -261.793,  -65.598,  265.354,  -2.131,  -265.039,
    -265.039,  -65.031,  261.039,  -2.128,  -260.869,
    -265.869,  -66.309,  265.897,  -2.727,  -265.958,
    -267.231,  -66.862,  267.255,  -2.545,  -267.009,
    -265.8  ,  -65.551,  266.254,  -2.394,  -264.825,
    -264.825,  -65.205,  262.245,  -2.312,  -268.365,
    -264.253,  -64.514,  264.413,  -2.413,  -264.413,
    -268.413,  -64.413,  264.589,  -2.277,  -265.378,
    -265.69 ,  -65.69 ,  264.972,  -2.972,  -269.972,
    -264.972,  -65.283,  265.237,  -2.671,  -264.88 ,
    -265.099,  -66.919,  265.878,  -2.653,  -264.653;

  Eigen::VectorXd all_draws(num_draws * num_params * num_chains);
  for (int i = 0; i < all_draws.size(); ++i) {
    all_draws(i) = draws_pool(i);
  }

  std::vector<Eigen::VectorXd> chain_draws(num_draws);
  for (int i = 0; i < num_draws; ++i) {
    int i_begin = comm.rank() * num_draws * num_params + i * num_params;
    chain_draws[i] = Eigen::VectorXd::Map(&all_draws(i_begin), num_params);
  }

  stan::math::mpi::mpi_covar_estimator estimator(num_params, num_draws);
  for (int i = 0; i < num_draws; ++i) {
    estimator.add_sample(chain_draws[i]);
  }

  stan::math::welford_covar_estimator welford(num_params);
  for (int i = 0; i < num_draws * num_chains; ++i) {
    welford.add_sample(Eigen::VectorXd::Map(&all_draws(i * num_params), num_params));
  }

  EXPECT_EQ(num_draws * comm.size(), estimator.num_samples(comm));

  Eigen::MatrixXd cv, cv0;
  estimator.sample_covariance(cv, comm);
  welford.sample_covariance(cv0);
  for (int i = 0; i < cv.size(); ++i) {
    EXPECT_FLOAT_EQ(cv(i), cv0(i));
  }
}

TEST(mpi_covariance_test, covar_use_recent_draws) {
  const Communicator comm(MPI_COMM_STAN);

  const int num_draws = 3;
  const int num_params = 5;
  const int max_num_chains = 4;
  const int num_chains = comm.size();
  Eigen::VectorXd draws_pool(num_draws * num_params * max_num_chains);
  draws_pool <<
    -266.132,  -65.76 ,  265.757,  -2.793,  -265.791,
    -261.793,  -65.598,  265.354,  -2.131,  -265.039,
    -265.039,  -65.031,  261.039,  -2.128,  -260.869,
    -265.869,  -66.309,  265.897,  -2.727,  -265.958,
    -267.231,  -66.862,  267.255,  -2.545,  -267.009,
    -265.8  ,  -65.551,  266.254,  -2.394,  -264.825,
    -264.825,  -65.205,  262.245,  -2.312,  -268.365,
    -264.253,  -64.514,  264.413,  -2.413,  -264.413,
    -268.413,  -64.413,  264.589,  -2.277,  -265.378,
    -265.69 ,  -65.69 ,  264.972,  -2.972,  -269.972,
    -264.972,  -65.283,  265.237,  -2.671,  -264.88 ,
    -265.099,  -66.919,  265.878,  -2.653,  -264.653;

  Eigen::VectorXd all_draws(num_draws * num_params * num_chains);
  for (int i = 0; i < all_draws.size(); ++i) {
    all_draws(i) = draws_pool(i);
  }

  std::vector<Eigen::VectorXd> chain_draws(num_draws);
  for (int i = 0; i < num_draws; ++i) {
    int i_begin = comm.rank() * num_draws * num_params + i * num_params;
    chain_draws[i] = Eigen::VectorXd::Map(&all_draws(i_begin), num_params);
  }

  stan::math::mpi::mpi_covar_estimator estimator(num_params, num_draws);
  for (int i = 0; i < num_draws; ++i) {
    estimator.add_sample(chain_draws[i]);
  }

  int col_begin = 1;
  int n_draw = num_draws - col_begin;;
  stan::math::welford_covar_estimator welford(num_params);
  for (int chain = 0; chain < num_chains; ++chain) {
    for (int i = 0; i < n_draw; ++i) {
      int j = chain * num_draws + col_begin + i;
      welford.add_sample(Eigen::VectorXd::Map(&all_draws(j * num_params), num_params));
    }
  }

  Eigen::MatrixXd cv, cv0;
  estimator.sample_covariance(cv, col_begin, n_draw, comm);
  welford.sample_covariance(cv0);
  for (int i = 0; i < cv.size(); ++i) {
    EXPECT_FLOAT_EQ(cv(i), cv0(i));
  }
}

#endif
