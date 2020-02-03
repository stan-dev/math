// the tests here check that map_rect_concurrent works correct as such we
// enforce that STAN_MPI is NOT defined

#ifdef STAN_MPI
#undef STAN_MPI
#endif

#include <gtest/gtest.h>
#include <stan/math/rev.hpp>

#include <test/unit/math/prim/functor/hard_work.hpp>
#include <test/unit/math/prim/functor/utils_threads.hpp>

#include <iostream>
#include <vector>

STAN_REGISTER_MAP_RECT(0, hard_work)
STAN_REGISTER_MAP_RECT(1, hard_work)

struct map_rect : public ::testing::Test {
  Eigen::VectorXd shared_params_d;
  std::vector<Eigen::VectorXd> job_params_d;
  std::vector<std::vector<double> > x_r;
  std::vector<std::vector<int> > x_i;
  const int N = 100;

  virtual void SetUp() {
    set_n_threads(4);
    shared_params_d.resize(2);
    shared_params_d << 2, 0;

    for (int n = 0; n != N; ++n) {
      Eigen::VectorXd job_d(2);
      job_d << 0, n * n;
      job_params_d.push_back(job_d);
    }

    x_r = std::vector<std::vector<double> >(N, std::vector<double>(1, 1.0));
    x_i = std::vector<std::vector<int> >(N, std::vector<int>(1, 0));
  }
};

TEST_F(map_rect, concurrent_ragged_return_size_dd) {
  Eigen::VectorXd res1 = stan::math::map_rect<0, hard_work>(
      shared_params_d, job_params_d, x_r, x_i);

  EXPECT_EQ(res1.size(), 2 * N);

  x_i[1] = std::vector<int>(1, 1);

  Eigen::VectorXd res2 = stan::math::map_rect<1, hard_work>(
      shared_params_d, job_params_d, x_r, x_i);

  EXPECT_EQ(res2.size(), 2 * N + 1);
}

TEST_F(map_rect, concurrent_eval_ok_dd) {
  Eigen::VectorXd res1 = stan::math::map_rect<0, hard_work>(
      shared_params_d, job_params_d, x_r, x_i);
  for (int i = 0, j = 0; i < N; i++) {
    j = 2 * i;
    EXPECT_FLOAT_EQ(
        res1(j), job_params_d[i](0) * job_params_d[i](0) + shared_params_d(0));
    EXPECT_FLOAT_EQ(res1(j + 1),
                    x_r[i][0] * job_params_d[i](1) * job_params_d[i](0)
                        + 2 * shared_params_d(0) + shared_params_d(1));
  }
}
