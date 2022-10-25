#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>
#include <vector>

namespace exp_mod_normal_lccdf4_test {

auto exp_mod_normal_lccdf_functor
    = [](const auto& y, const auto& mu, const auto& sigma, const auto& lambda) {
        return stan::math::exp_mod_normal_lccdf(y, mu, sigma, lambda);
      };

TEST(ProbDistributionsDoubleExpModNormalLccdf, opencl_broadcast_lambda) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0.3, 0.8, 1.0;
  Eigen::VectorXd mu(N);
  mu << 0.3, 0.8, 1.0;
  Eigen::VectorXd sigma(N);
  sigma << 0.3, 0.4, 1.1;
  double lambda_scal = 12.3;

  stan::math::test::test_opencl_broadcasting_prim_rev<3>(
      exp_mod_normal_lccdf_functor, y, mu, sigma, lambda_scal);
  stan::math::test::test_opencl_broadcasting_prim_rev<3>(
      exp_mod_normal_lccdf_functor, y.transpose().eval(), mu,
      sigma.transpose().eval(), lambda_scal);
}
}  // namespace exp_mod_normal_lccdf4_test

#endif
