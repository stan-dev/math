#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>
#include <vector>

auto exp_mod_normal_lcdf_functor
    = [](const auto& y, const auto& mu, const auto& sigma, const auto& lambda) {
        return stan::math::exp_mod_normal_lcdf(y, mu, sigma, lambda);
      };

TEST(ProbDistributionsDoubleExpModNormalLcdf, opencl_broadcast_mu) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0.3, 0.8, 1.0;
  double mu_scal = 12.3;
  Eigen::VectorXd sigma(N);
  sigma << 0.3, 0.8, 1.0;
  Eigen::VectorXd lambda(N);
  lambda << 0.3, 0.4, 1.1;

  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      exp_mod_normal_lcdf_functor, y, mu_scal, sigma, lambda);
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      exp_mod_normal_lcdf_functor, y.transpose().eval(), mu_scal, sigma,
      lambda.transpose().eval());
}

#endif
