#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>
#include <gtest/gtest.h>

TEST_F(AgradRev, logistic_lccdf_upper_tail) {
  for (double y_value : {36.5, 40.0}) {
    stan::math::var y = y_value;
    stan::math::var mu = 0.0;
    stan::math::var sigma = 1.0;

    stan::math::var log_ccdf = stan::math::logistic_lccdf(y, mu, sigma);
    log_ccdf.grad();

    const double deriv = stan::math::inv_logit(y_value);
    EXPECT_DOUBLE_EQ(-stan::math::log1p_exp(y_value), log_ccdf.val());
    EXPECT_DOUBLE_EQ(-deriv, y.adj());
    EXPECT_DOUBLE_EQ(deriv, mu.adj());
    EXPECT_DOUBLE_EQ(y_value * deriv, sigma.adj());
    stan::math::recover_memory();
  }
}
