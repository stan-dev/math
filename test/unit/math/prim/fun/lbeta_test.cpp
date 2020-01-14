#include <test/unit/math/expect_near_rel.hpp>
#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <vector>

TEST(MathFunctions, lbeta) {
  using stan::math::lbeta;

  EXPECT_FLOAT_EQ(0.0, lbeta(1.0, 1.0));
  EXPECT_FLOAT_EQ(2.981361, lbeta(0.1, 0.1));
  EXPECT_FLOAT_EQ(-4.094345, lbeta(3.0, 4.0));
  EXPECT_FLOAT_EQ(-4.094345, lbeta(4.0, 3.0));
}

TEST(MathFunctions, lbeta_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_TRUE(std::isnan(stan::math::lbeta(nan, 1.0)));

  EXPECT_TRUE(std::isnan(stan::math::lbeta(1.0, nan)));

  EXPECT_TRUE(std::isnan(stan::math::lbeta(nan, nan)));
}

TEST(MathFunctions, lbeta_extremes_errors) {
  double inf = std::numeric_limits<double>::infinity();
  double after_stirling = 
    std::nextafter(stan::math::lgamma_stirling_diff_useful, inf);
  using stan::math::lbeta;

  EXPECT_FLOAT_EQ(lbeta(0.0, 1.0), inf);
  EXPECT_FLOAT_EQ(lbeta(1.0, 0.0), inf);
  EXPECT_FLOAT_EQ(lbeta(0.0, after_stirling), inf);
  EXPECT_FLOAT_EQ(lbeta(after_stirling, 0.0), inf);
  EXPECT_FLOAT_EQ(lbeta(0.0,0.0), inf);

  EXPECT_FLOAT_EQ(lbeta(inf, 0.0), inf);
  EXPECT_FLOAT_EQ(lbeta(0.0, inf), inf);
  EXPECT_FLOAT_EQ(lbeta(inf, 1), -inf);
  EXPECT_FLOAT_EQ(lbeta(1e8, inf), -inf);
  EXPECT_FLOAT_EQ(lbeta(inf, inf), -inf);
}

TEST(MathFunctions, lbeta_identities) {
  using stan::test::expect_near_rel;
  using stan::math::lbeta;
  using stan::math::pi;

  std::vector<double> to_test = 
    {1e-100, 1e-8, 1e-1, 1, 1 + 1e-6, 1e3, 1e30, 1e100};
  for(double x : to_test) {
    for(double y : to_test) {
      std::ostringstream msg;
      msg << "successors: x = " << x << "; y = " << y;
      expect_near_rel(msg.str(), lbeta(x, y), 
        stan::math::log_sum_exp(lbeta(x + 1, y), lbeta(x, y + 1)));
    }
  }

  for(double x : to_test) {
    if(x < 1) {
      std::ostringstream msg;
      msg << "sin: x = " << x ;
      expect_near_rel(msg.str(), lbeta(x, 1.0 - x), 
        log(pi()) - log(sin(pi() * x)));
    }
  }

  for(double x : to_test) {
    std::ostringstream msg;
    msg << "inv: x = " << x ;
    expect_near_rel(msg.str(), lbeta(x, 1.0), 
      -log(x));
  }
}

TEST(MathFunctions, lbeta_stirling_cutoff) {
  using stan::test::expect_near_rel;

  double after_stirling = 
    std::nextafter(stan::math::lgamma_stirling_diff_useful, 
      std::numeric_limits<double>::infinity());
  double before_stirling = 
    std::nextafter(stan::math::lgamma_stirling_diff_useful, 0);
  using stan::math::lbeta;

  std::vector<double> to_test = 
    {1e-100, 1e-8, 1e-1, 1, 1 + 1e-6, 1e3, 1e30, 1e100, 
    before_stirling, after_stirling};
  for(double x : to_test) {
    std::ostringstream msg;
    msg << "before and after cutoff: x = " << x 
      << "; cutoff = " << stan::math::lgamma_stirling_diff_useful;
    expect_near_rel(msg.str(), 
      lbeta(x, before_stirling), lbeta(x, after_stirling));
    expect_near_rel(msg.str(), 
      lbeta(before_stirling, x), lbeta(after_stirling, x));
  }
}