#ifdef STAN_OPENCL
#ifndef STAN_TEST_SKIP_REQUIRING_OPENCL_INT64_BASE_ATOMIC
#include <stan/math/opencl/indexing_rev.hpp>
#include <stan/math/opencl/copy.hpp>
#include <test/unit/util.hpp>
#include <test/unit/math/expect_near_rel.hpp>
#include <gtest/gtest.h>
#include <algorithm>

TEST(indexing_rev, indexing_rev_small) {
  Eigen::MatrixXd res(3, 3);
  res << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  Eigen::MatrixXi idx(3, 3);
  idx << 1, 2, 1, 0, 1, 2, 1, 1, 1;
  Eigen::MatrixXd adj = Eigen::MatrixXd::Zero(3, 1);
  Eigen::MatrixXd correct(3, 1);
  correct << 4, 33, 8;
  stan::math::matrix_cl<double> res_cl(res);
  stan::math::matrix_cl<int> idx_cl(idx);
  stan::math::matrix_cl<double> adj_cl(adj);

  stan::math::indexing_rev(adj_cl, idx_cl, res_cl);
  EXPECT_NEAR_REL(stan::math::from_matrix_cl(adj_cl), correct);
}

TEST(indexing_rev, indexing_rev_large) {
  int N = 377;

  // different sizes of adj ensure all three kernels are tested regardless of
  // the OpenCL device
  for (int M = 1; M < 1e6; M *= 2) {
    Eigen::MatrixXd res = Eigen::MatrixXd::Random(N, N);
    Eigen::MatrixXi idx(N, N);
    for (int i = 0; i < N * N; i++) {
      idx(i) = std::abs(Eigen::MatrixXi::Random(1, 1)(0)) % M;
    }
    Eigen::MatrixXd adj = Eigen::MatrixXd::Zero(M, 1);
    Eigen::MatrixXd correct = adj;
    for (int i = 0; i < N * N; i++) {
      correct(idx(i)) += res(i);
    }
    stan::math::matrix_cl<double> res_cl(res);
    stan::math::matrix_cl<int> idx_cl(idx);
    stan::math::matrix_cl<double> adj_cl(adj);

    stan::math::indexing_rev(adj_cl, idx_cl, res_cl);
    EXPECT_NEAR_REL(stan::math::from_matrix_cl(adj_cl), correct);
  }
}
#endif  // STAN_TEST_SKIP_REQUIRING_OPENCL_INT64_BASE_ATOMIC
#endif
