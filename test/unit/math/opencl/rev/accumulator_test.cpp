#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>

TEST(OpenCLAcc, var_mat) {
  auto f = [](const auto& x) {
    using x_t = decltype(x);
    stan::math::accumulator<stan::scalar_type_t<x_t>> acc;
    acc.add(x);
    return acc.sum();
  };
  Eigen::MatrixXd x = Eigen::MatrixXd::Random(2, 2);
  stan::math::test::compare_cpu_opencl_prim_rev(f, x);
  std::vector<Eigen::MatrixXd> x_vec;
  x_vec.push_back(x);
  x_vec.push_back(x);
  x_vec.push_back(x);
  stan::math::test::compare_cpu_opencl_prim_rev(f, x_vec);
}
#endif
