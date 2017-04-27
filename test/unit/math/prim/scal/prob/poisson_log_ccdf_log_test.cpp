#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

void check_values(int count, double log_rate) {
  using std::exp;
  using stan::math::is_nan;
  using stan::math::poisson_log_ccdf_log;
  using stan::math::poisson_ccdf_log;

  double rate = exp(log_rate);

  double val1 = poisson_log_ccdf_log(count, log_rate);
  double val2 = poisson_ccdf_log(count, rate);

  EXPECT_FALSE(is_nan(val1));
  EXPECT_FALSE(val1 > 0.0);
  EXPECT_FLOAT_EQ(val1, val2);
}

TEST(ProbPoisson, log_ccdf_log_values) {
  int count = 2;
  for (int i = 0; i < 6; ++i) {
    double log_rate = 0.3;
    for (int j = 0; j < 5; ++j) {
      check_values(count, log_rate);
    }
  }
}
