#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>

#include <test/unit/math/prim/mat/functor/hard_work.hpp>

#include <iostream>
#include <vector>

// the tests here check that map_rect refuses mal-formatted input as
// such it does not matter if STAN_MPI is defined or not

STAN_REGISTER_MAP_RECT(0, hard_work)
STAN_REGISTER_MAP_RECT(1, hard_work)
STAN_REGISTER_MAP_RECT(2, hard_work)
STAN_REGISTER_MAP_RECT(3, hard_work)
STAN_REGISTER_MAP_RECT(4, hard_work)
STAN_REGISTER_MAP_RECT(5, hard_work)

struct map_rect : public ::testing::Test {
  Eigen::VectorXd shared_params_d;
  std::vector<Eigen::VectorXd> job_params_d;
  std::vector<std::vector<double> > x_r;
  std::vector<std::vector<int> > x_i;
  const std::size_t N = 10;

  virtual void SetUp() {
    shared_params_d.resize(2);
    shared_params_d << 2, 0;

    for (std::size_t n = 0; n != N; ++n) {
      Eigen::VectorXd job_d(2);
      job_d << 0, n * n;
      job_params_d.push_back(job_d);
    }

    x_r = std::vector<std::vector<double> >(N, std::vector<double>(1, 1.0));
    x_i = std::vector<std::vector<int> >(N, std::vector<int>(1, 0));
  }
};

TEST_F(map_rect, no_job_input_ok_dd) {
  job_params_d.resize(0);
  x_i.resize(0);
  x_r.resize(0);

  EXPECT_NO_THROW((stan::math::map_rect<0, hard_work>(shared_params_d,
                                                      job_params_d, x_r, x_i)));
}

TEST_F(map_rect, size_mismatch_job_params_dd) {
  job_params_d.pop_back();

  EXPECT_THROW((stan::math::map_rect<1, hard_work>(shared_params_d,
                                                   job_params_d, x_r, x_i)),
               std::invalid_argument);
}

TEST_F(map_rect, size_mismatch_real_data_dd) {
  x_r.pop_back();

  EXPECT_THROW((stan::math::map_rect<2, hard_work>(shared_params_d,
                                                   job_params_d, x_r, x_i)),
               std::invalid_argument);
}

TEST_F(map_rect, size_mismatch_int_data_dd) {
  x_i.pop_back();

  EXPECT_THROW((stan::math::map_rect<3, hard_work>(shared_params_d,
                                                   job_params_d, x_r, x_i)),
               std::invalid_argument);
}

TEST_F(map_rect, wrong_size_job_params_dd) {
  job_params_d[1].resize(5);

  EXPECT_THROW((stan::math::map_rect<1, hard_work>(shared_params_d,
                                                   job_params_d, x_r, x_i)),
               std::invalid_argument);
}

TEST_F(map_rect, wrong_size_real_data_dd) {
  x_r[1] = std::vector<double>(5, 1);

  EXPECT_THROW((stan::math::map_rect<4, hard_work>(shared_params_d,
                                                   job_params_d, x_r, x_i)),
               std::invalid_argument);
}

TEST_F(map_rect, wrong_size_int_data_dd) {
  x_i[1] = std::vector<int>(5, 1);

  EXPECT_THROW((stan::math::map_rect<5, hard_work>(shared_params_d,
                                                   job_params_d, x_r, x_i)),
               std::invalid_argument);
}
