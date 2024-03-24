#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>

TEST(mathMixScalFun, weibull_lcdf) {
  // Inputs are tested on the log (i.e., unconstrained) scale so that the
  // finite-diffs don't result in invalid inputs.
  auto f = [](const auto& y, const auto& alpha, const auto& sigma) {
    using stan::math::lb_constrain;
    using stan::math::positive_constrain;

    return stan::math::weibull_lcdf(lb_constrain(y, 0.0),
                                    positive_constrain(alpha),
                                    positive_constrain(sigma));
  };

  using stan::math::log;

  Eigen::VectorXd y(3);
  y << stan::math::NEGATIVE_INFTY, 1.2, 0.0;  // lb_constrain(y[0], 0.0) = 0.0

  Eigen::VectorXd alpha(3);
  alpha << 2.0, 3.0, 4.0;

  Eigen::VectorXd sigma(3);
  sigma << 5.0, 6.0, 7.0;

  stan::test::expect_ad(f, y, alpha, sigma);
  stan::test::expect_ad(f, y[0], alpha, sigma);
  stan::test::expect_ad(f, y, alpha[0], sigma);
  stan::test::expect_ad(f, y, alpha, sigma[0]);

  stan::test::expect_ad(f, y[0], alpha[0], sigma);
  stan::test::expect_ad(f, y[0], alpha, sigma[0]);
  stan::test::expect_ad(f, y, alpha[0], sigma[0]);

  stan::test::expect_ad(f, y[0], alpha[0], sigma[0]);
}
