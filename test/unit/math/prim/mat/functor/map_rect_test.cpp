#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>

#include <test/unit/math/prim/mat/functor/hard_work.hpp>

#include <iostream>

TEST(map_rect_serial, hard_work_dd) {
  Eigen::VectorXd shared_params_d(2);
  shared_params_d << 2, 0;
  std::vector<Eigen::VectorXd> job_params_d;

  const std::size_t N = 10;

  for(std::size_t n = 0; n != N; ++n) {
    Eigen::VectorXd job_d(2);
    job_d << 0, n * n;
    job_params_d.push_back(job_d);
  }

  std::vector<std::vector<double> > x_r(N, std::vector<double>(1,1.0));
  std::vector<std::vector<int> > x_i(N, std::vector<int>(0));

  Eigen::VectorXd result = stan::math::map_rect<hard_work>(shared_params_d, job_params_d, x_r, x_i, 0);
}


