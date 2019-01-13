// the tests here check that map_rect_concurrent works correct as such we
// enforce that STAN_MPI is NOT defined and these tests only run if
// STAN_THREADS is defined

#ifdef STAN_THREADS

#ifdef STAN_MPI
#undef STAN_MPI
#endif

#include <stdlib.h>

#include <gtest/gtest.h>
#include <stan/math/prim/mat.hpp>

#include <test/unit/math/prim/mat/functor/hard_work.hpp>

#include <iostream>
#include <vector>
#include <string>

// utility to set number of threads to use
void set_n_threads(int num_threads) {
  static char env_string[256];
  std::string num_threads_str = std::to_string(num_threads);
  snprintf(env_string, sizeof(env_string), "STAN_NUM_THREADS=%s",
           num_threads_str.c_str());
  putenv(env_string);
}

STAN_REGISTER_MAP_RECT(0, hard_work)
STAN_REGISTER_MAP_RECT(1, hard_work)

void setup_job(int N, Eigen::VectorXd& shared_params_d,
               std::vector<Eigen::VectorXd>& job_params_d,
               std::vector<std::vector<double> >& x_r,
               std::vector<std::vector<int> >& x_i) {
  shared_params_d.resize(2);
  shared_params_d << 2, 0;

  job_params_d.clear();
  for (int n = 0; n != N; ++n) {
    Eigen::VectorXd job_d(2);
    job_d << 1.1, n * n;
    job_params_d.push_back(job_d);
  }

  x_r.clear();
  x_i.clear();

  if (N > 0) {
    x_r = std::vector<std::vector<double> >(N, std::vector<double>(1, 1.0));
    x_i = std::vector<std::vector<int> >(N, std::vector<int>(1, 0));
  }
}

struct map_rect : public ::testing::Test {
  Eigen::VectorXd shared_params_d;
  std::vector<Eigen::VectorXd> job_params_d;
  std::vector<std::vector<double> > x_r;
  std::vector<std::vector<int> > x_i;
  std::vector<int> N_test{0, 1, 2, 3, 4, 5, 7, 10, 13, 67, 100};

  virtual void SetUp() {}
};

TEST_F(map_rect, concurrent_varying_num_threads_ragged_dd) {
  set_n_threads(-1);

  for (std::size_t i = 1; i < 20; ++i) {
    for (std::size_t n = 0; n < N_test.size(); ++n) {
      const int N = N_test[n];

      setup_job(N, shared_params_d, job_params_d, x_r, x_i);

      Eigen::VectorXd res1 = stan::math::map_rect<0, hard_work>(
          shared_params_d, job_params_d, x_r, x_i);
      EXPECT_EQ(res1.size(), 2 * N);
      for (int i = 0, j = 0; i < N; i++) {
        j = 2 * i;
        EXPECT_FLOAT_EQ(res1(j), job_params_d[i](0) * job_params_d[i](0)
                                     + shared_params_d(0));
        EXPECT_FLOAT_EQ(res1(j + 1),
                        x_r[i][0] * job_params_d[i](1) * job_params_d[i](0)
                            + 2 * shared_params_d(0) + shared_params_d(1));
      }
    }
    set_n_threads(i);
  }

  set_n_threads(-1);
  for (std::size_t i = 1; i < 20; ++i) {
    for (std::size_t n = 1; n < N_test.size(); ++n) {
      const int N = N_test[n];
      setup_job(N, shared_params_d, job_params_d, x_r, x_i);

      x_i[0] = std::vector<int>(1, 1);

      Eigen::VectorXd res2 = stan::math::map_rect<0, hard_work>(
          shared_params_d, job_params_d, x_r, x_i);
      EXPECT_EQ(res2.size(), 2 * N + 1);
    }
    set_n_threads(i);
  }
}

TEST_F(map_rect, concurrent_varying_num_threads_eval_ok_dd) {
  set_n_threads(-1);
  for (std::size_t i = 1; i < 20; ++i) {
    for (std::size_t n = 0; n < N_test.size(); ++n) {
      const int N = N_test[n];
      setup_job(N, shared_params_d, job_params_d, x_r, x_i);

      Eigen::VectorXd res1 = stan::math::map_rect<0, hard_work>(
          shared_params_d, job_params_d, x_r, x_i);
      EXPECT_EQ(res1.size(), 2 * N);
      for (int i = 0, j = 0; i < N; i++) {
        j = 2 * i;
        EXPECT_FLOAT_EQ(res1(j), job_params_d[i](0) * job_params_d[i](0)
                                     + shared_params_d(0));
        EXPECT_FLOAT_EQ(res1(j + 1),
                        x_r[i][0] * job_params_d[i](1) * job_params_d[i](0)
                            + 2 * shared_params_d(0) + shared_params_d(1));
      }
    }
    set_n_threads(i);
  }
}

#endif
