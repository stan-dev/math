#include <stan/math/prim.hpp>
#include <test/unit/math/prim/fun/binary_scalar_tester.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <string>
#include <vector>
#include <algorithm>

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
  double after_stirling
      = std::nextafter(stan::math::lgamma_stirling_diff_useful, inf);
  using stan::math::lbeta;

  EXPECT_FLOAT_EQ(lbeta(0.0, 1.0), inf);
  EXPECT_FLOAT_EQ(lbeta(1.0, 0.0), inf);
  EXPECT_FLOAT_EQ(lbeta(0.0, after_stirling), inf);
  EXPECT_FLOAT_EQ(lbeta(after_stirling, 0.0), inf);
  EXPECT_FLOAT_EQ(lbeta(0.0, 0.0), inf);

  EXPECT_FLOAT_EQ(lbeta(inf, 0.0), inf);
  EXPECT_FLOAT_EQ(lbeta(0.0, inf), inf);
  EXPECT_FLOAT_EQ(lbeta(inf, 1), -inf);
  EXPECT_FLOAT_EQ(lbeta(1e8, inf), -inf);
  EXPECT_FLOAT_EQ(lbeta(inf, inf), -inf);
}

TEST(MathFunctions, lbeta_identities) {
  using stan::math::lbeta;
  using stan::math::pi;

  std::vector<double> to_test
      = {1e-100, 1e-8, 1e-1, 1, 1 + 1e-6, 1e3, 1e30, 1e100};
  auto tol = [](double x, double y) {
    return std::max(1e-15 * (0.5 * (fabs(x) + fabs(y))), 1e-15);
  };

  for (double x : to_test) {
    for (double y : to_test) {
      std::stringstream msg;
      msg << std::setprecision(22) << "successors: x = " << x << "; y = " << y;
      double lh = lbeta(x, y);
      double rh = stan::math::log_sum_exp(lbeta(x + 1, y), lbeta(x, y + 1));
      EXPECT_NEAR(lh, rh, tol(lh, rh)) << msg.str();
    }
  }

  for (double x : to_test) {
    if (x < 1) {
      std::stringstream msg;
      msg << std::setprecision(22) << "sin: x = " << x;
      double lh = lbeta(x, 1.0 - x);
      double rh = log(pi()) - log(sin(pi() * x));
      EXPECT_NEAR(lh, rh, tol(lh, rh)) << msg.str();
    }
  }

  for (double x : to_test) {
    std::stringstream msg;
    msg << std::setprecision(22) << "inv: x = " << x;
    double lh = lbeta(x, 1.0);
    double rh = -log(x);
    EXPECT_NEAR(lh, rh, tol(lh, rh)) << msg.str();
  }
}

TEST(MathFunctions, lbeta_stirling_cutoff) {
  using stan::math::lgamma_stirling_diff_useful;

  double after_stirling
      = std::nextafter(lgamma_stirling_diff_useful, stan::math::INFTY);
  double before_stirling = std::nextafter(lgamma_stirling_diff_useful, 0);
  using stan::math::lbeta;

  std::vector<double> to_test
      = {1e-100,          1e-8,          1e-1, 1, 1 + 1e-6, 1e3, 1e30, 1e100,
         before_stirling, after_stirling};
  for (const double x : to_test) {
    double before = lbeta(x, before_stirling);
    double at = lbeta(x, lgamma_stirling_diff_useful);
    double after = lbeta(x, after_stirling);

    double diff_before = at - before;
    double diff_after = after - at;
    double tol
        = std::max(1e-15 * (0.5 * (fabs(diff_before) + fabs(diff_after))),
                   1e-14 * fabs(at));

    EXPECT_NEAR(diff_before, diff_after, tol)
        << "diff before and after cutoff: x = " << x << "; before = " << before
        << "; at = " << at << "; after = " << after;
  }
}

TEST(MathFunctions, lbeta_vec) {
  auto f = [](const auto& x1, const auto& x2) {
    using stan::math::lbeta;
    return lbeta(x1, x2);
  };

  Eigen::VectorXd in1(3);
  in1 << 1.8, 3.24, 1.8;
  Eigen::VectorXd in2(3);
  in2 << 1.3, 0.7, 2.8;
  stan::test::binary_scalar_tester(f, in1, in2);
}
