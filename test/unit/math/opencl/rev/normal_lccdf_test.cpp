#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>
#include <vector>


auto normal_lccdf_functor
    = [](const auto& y, const auto& mu, const auto& sigma) {
        return stan::math::normal_lccdf(y, mu, sigma);
      };



TEST(ProbDistributionsNormalLccdf, opencl_matches_cpu_big) {
  int N = 153;

std::srand(123);
for (int i = 0; i < 10; ++i) {
  Eigen::Matrix<double, Eigen::Dynamic, 1> mu
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1) + 1.0;
  Eigen::Matrix<double, Eigen::Dynamic, 1> sigma
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs() + 0.01;
  Eigen::Matrix<double, Eigen::Dynamic, 1> y = (mu.array() * sigma.array()).matrix();
  std::cout << "Iter: " << i << " mu, sigma, y" << std::endl;
  for (int j = 0; j < N; j++) {
    std::cout << mu(j) << ", " << sigma(j) << ", " << y(j) << std::endl;
  }
  std::cout << "-----------compare_cpu_opencl_prim_rev" << std::endl;
  stan::math::test::compare_cpu_opencl_prim_rev(normal_lccdf_functor, y, mu,
                                                sigma);
  std::cout << "-----------compare_cpu_opencl_prim_rev transpose" << std::endl;
  stan::math::test::compare_cpu_opencl_prim_rev(
      normal_lccdf_functor, y.transpose().eval(), mu.transpose().eval(),
      sigma.transpose().eval());
}
}
#endif
