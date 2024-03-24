#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>

TEST(mathMixScalFun, weibull_cdf) {
  // Inputs are tested on the log (i.e., unconstrained) scale so that the
  // finite-diffs don't result in invalid inputs.
  auto f = [](const auto& y, const auto& alpha, const auto& sigma) {
    using stan::math::lb_constrain;
    using stan::math::positive_constrain;

    return stan::math::weibull_cdf(lb_constrain(y, 0.0),
                                    positive_constrain(alpha),
                                    positive_constrain(sigma));
  };

  using stan::math::log;

  stan::test::expect_ad(f, log(1.2), log(4), log(20));
  stan::test::expect_ad(f, 0.0, 0.0, 0.0);
  stan::test::expect_ad(f, stan::math::NEGATIVE_INFTY, 0.0, 0.0); // y = 0.0
}
