#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(ProbPoisson, log_cdf_log_chain_rule) {
  using std::exp;

  using stan::math::is_nan;
  using stan::math::var;

  using stan::math::poisson_log_cdf_log;
  using stan::math::poisson_cdf_log;

  int count = 20;
  double log_rate_dbl = 2.1;
  double rate_dbl = exp(log_rate_dbl);
  
  var log_rate(log_rate_dbl);
  var rate(rate_dbl);

  var val1 = poisson_log_cdf_log(count, log_rate);
  var val2 = poisson_cdf_log(count, rate);
  
  std::vector<var> x1;
  x1.push_back(log_rate);

  std::vector<var> x2;
  x2.push_back(rate);

  std::vector<double> gradients1;
  std::vector<double> gradients2;

  val1.grad(x1, gradients1);
  val2.grad(x2, gradients2);

  EXPECT_FALSE(is_nan(gradients1[0]));
  EXPECT_FALSE(is_nan(gradients2[0]));
  EXPECT_NEAR(gradients1[0], gradients2[0] * rate_dbl, 0.01);
}
