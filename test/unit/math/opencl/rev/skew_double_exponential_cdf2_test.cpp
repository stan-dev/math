#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>
#include <vector>

namespace skew_double_exponential_cdf2_test {

auto skew_double_exponential_cdf_functor
    = [](const auto& y, const auto& mu, const auto& sigma, const auto& tau) {
        return stan::math::skew_double_exponential_cdf(y, mu, sigma, tau);
      };

TEST(ProbDistributionsSkewDoubleExponentialCdf, opencl_broadcast_sigma) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0.3, 0.8, 1.0;
  Eigen::VectorXd mu(N);
  mu << 0.3, 0.8, 1.0;
  double sigma_scal = 12.3;
  Eigen::VectorXd tau(N);
  tau << 0.3, 0.4, 0.9;

  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      skew_double_exponential_cdf_functor, y, mu, sigma_scal, tau);
  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      skew_double_exponential_cdf_functor, y.transpose().eval(), mu, sigma_scal,
      tau.transpose().eval());
}

TEST(ProbDistributionsSkewDoubleExponentialCdf, opencl_broadcast_tau) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0.3, 0.8, 1.0;
  Eigen::VectorXd mu(N);
  mu << 0.3, 0.8, 1.0;
  Eigen::VectorXd sigma(N);
  sigma << 0.3, 0.4, 1.1;
  double tau_scal = 0.3;

  stan::math::test::test_opencl_broadcasting_prim_rev<3>(
      skew_double_exponential_cdf_functor, y, mu, sigma, tau_scal);
  stan::math::test::test_opencl_broadcasting_prim_rev<3>(
      skew_double_exponential_cdf_functor, y.transpose().eval(), mu,
      sigma.transpose().eval(), tau_scal);
}
}  // namespace skew_double_exponential_cdf2_test

#endif
