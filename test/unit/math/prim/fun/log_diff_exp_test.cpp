#include <stan/math/prim.hpp>
#include <test/unit/math/prim/fun/binary_scalar_tester.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

void test_log_diff_exp(double a, double b) {
  using stan::math::log_diff_exp;
  using std::exp;
  using std::log;
  EXPECT_FLOAT_EQ(log(exp(a) - exp(b)), log_diff_exp(a, b));
}

TEST(MathFunctions, log_diff_exp) {
  using stan::math::log_diff_exp;
  test_log_diff_exp(3.0, 2.0);
  test_log_diff_exp(4.0, 1.0);
  test_log_diff_exp(3.0, 2.0);
  test_log_diff_exp(0, -2.1);
  test_log_diff_exp(-20.0, -23);
  test_log_diff_exp(-21.2, -32.1);
  EXPECT_NO_THROW(log_diff_exp(-20.0, 12));
  EXPECT_NO_THROW(log_diff_exp(-20.0, -12.1));
  EXPECT_NO_THROW(log_diff_exp(120.0, 120.10));
  EXPECT_NO_THROW(log_diff_exp(-20.0, 10.2));
  EXPECT_NO_THROW(log_diff_exp(10, 11));
  EXPECT_NO_THROW(log_diff_exp(10, 10));
  EXPECT_NO_THROW(log_diff_exp(-10.21, -10.21));

  // exp(10000.0) overflows
  EXPECT_FLOAT_EQ(10000.0, log_diff_exp(10000.0, 0.0));
  EXPECT_FLOAT_EQ(0.0, log_diff_exp(0.0, -10000.0));
}

TEST(MathFunctions, log_diff_exp_inf) {
  using stan::math::log_diff_exp;
  double inf = std::numeric_limits<double>::infinity();
  test_log_diff_exp(0.0, 0.0);
  test_log_diff_exp(-10.21, -10.21);
  test_log_diff_exp(0.0, -inf);
  test_log_diff_exp(-inf, -inf);
  EXPECT_FLOAT_EQ(inf, log_diff_exp(inf, 3.0));
  EXPECT_TRUE(std::isnan(log_diff_exp(inf, inf)));
}

TEST(MathFunctions, log_diff_exp_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_TRUE(std::isnan(stan::math::log_diff_exp(3.0, nan)));

  EXPECT_TRUE(std::isnan(stan::math::log_diff_exp(nan, 2.0)));

  EXPECT_TRUE(std::isnan(stan::math::log_diff_exp(nan, nan)));
}

TEST(MathFunctions, log_diff_exp_vec) {
  auto f = [](const auto& x1, const auto& x2) {
    using stan::math::log_diff_exp;
    return log_diff_exp(x1, x2);
  };

  Eigen::VectorXd in1(3);
  in1 << 4.1, 3.24, 6.8;
  Eigen::VectorXd in2(3);
  in2 << 2.8, 1.7, 3.1;
  stan::test::binary_scalar_tester(f, in1, in2);
}
