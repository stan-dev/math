#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>
#include <vector>

namespace exp_mod_normal_lcdf3_test {

auto exp_mod_normal_lcdf_functor
    = [](const auto& y, const auto& mu, const auto& sigma, const auto& lambda) {
        return stan::math::exp_mod_normal_lcdf(y, mu, sigma, lambda);
      };

TEST(ProbDistributionsDoubleExpModNormalLcdf, opencl_broadcast_sigma) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0.3, 0.8, 1.0;
  Eigen::VectorXd mu(N);
  mu << 0.3, 0.8, 1.0;
  double sigma_scal = 12.3;
  Eigen::VectorXd lambda(N);
  lambda << 0.3, 0.4, 1.1;

  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      exp_mod_normal_lcdf_functor, y, mu, sigma_scal, lambda);
  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      exp_mod_normal_lcdf_functor, y.transpose().eval(), mu, sigma_scal,
      lambda.transpose().eval());
}
}  // namespace exp_mod_normal_lcdf3_test

#endif
