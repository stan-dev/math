#include <stan/math/fwd.hpp>
#include <gtest/gtest.h>

#include <cmath>

TEST(AgradFwd, log_sum_exp_derivative_does_not_overflow) {
  using stan::math::fvar;
  using stan::math::log_sum_exp;
  
  Eigen::Matrix<fvar<double>, Eigen::Dynamic, 1> x(3);
  x(0) = fvar<double>(1000.0, 1.0);
  x(1) = fvar<double>(1001.0, 2.0);
  x(2) = fvar<double>(999.0, 3.0);
  
  auto y = log_sum_exp(x);
  // log_sum_exp([1000, 1001, 999])
  //
  // = 1001 + log(exp(-1) + 1 + exp(-2)).
  const double expected_value
      = 1001.0 + std::log(std::exp(-1.0) + 1.0 + std::exp(-2.0));
  
  EXPECT_NEAR(expected_value, y.val(), 1e-12);
  EXPECT_TRUE(std::isfinite(y.val()));
  
  // The derivative is
  //
  // softmax(x) dot x.d().
  //
  // The old implementation formed exp(x) directly in the derivative,
  // which gives inf / inf for these values.
  
  const double denom = std::exp(-1.0) + 1.0 + std::exp(-2.0);
  
  const double p0 = std::exp(-1.0) / denom;
  const double p1 = 1.0 / denom;
  const double p2 = std::exp(-2.0) / denom;
  
  const double expected_derivative = p0 + 2.0 * p1 + 3.0 * p2;
  
  EXPECT_TRUE(std::isfinite(y.d()));
  EXPECT_NEAR(expected_derivative, y.d(), 1e-12);
}
