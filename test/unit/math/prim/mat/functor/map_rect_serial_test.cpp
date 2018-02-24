// specific tests of the map_rect_serial implementation

#ifdef STAN_HAS_MPI
#undef STAN_HAS_MPI
#endif

#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>

#include <test/unit/math/prim/mat/functor/hard_work.hpp>

#include <iostream>
#include <vector>

// the tests here check that map_rect refuses mal-formatted input as
// such it does not matter if STAN_HAS_MPI is defined or not

STAN_REGISTER_MAP_RECT(0, hard_work)
STAN_REGISTER_MAP_RECT(1, hard_work)

struct map_rect : public ::testing::Test {
  Eigen::VectorXd shared_params_d;
  std::vector<Eigen::VectorXd> job_params_d;
  std::vector<std::vector<double> > x_r;
  std::vector<std::vector<int> > x_i;
  const int N = 10;

  virtual void SetUp() {
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


TEST_F(map_rect, ragged_return_size_dd) {
  Eigen::VectorXd res1 = stan::math::map_rect<0, hard_work>(shared_params_d,
                                                            job_params_d, x_r, x_i);

  EXPECT_EQ(res1.size(), 2*N);

  x_i[1] = std::vector<int>(1, 1);
  
  Eigen::VectorXd res2 = stan::math::map_rect<1, hard_work>(shared_params_d,
                                                            job_params_d, x_r, x_i);

  EXPECT_EQ(res2.size(), 2*N+1);
}
