#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

namespace {
// d/dx erfcx(x) = 2 * x * erfcx(x) - 2 / sqrt(pi)
double expected_deriv(double x) {
  return 2.0 * x * stan::math::erfcx(x) - stan::math::TWO_OVER_SQRT_PI;
}
}  // namespace

TEST_F(AgradRev, erfcx_value_and_gradient) {
  for (double x_value : {-2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0, 5.0}) {
    stan::math::var x = x_value;
    stan::math::var y = stan::math::erfcx(x);
    y.grad();
    EXPECT_FLOAT_EQ(stan::math::erfcx(x_value), y.val());
    EXPECT_FLOAT_EQ(expected_deriv(x_value), x.adj());
    stan::math::recover_memory();
  }
}

// Reference derivatives from 2 * x * erfcx(x) - 2 / sqrt(pi) evaluated in
// long double.
TEST_F(AgradRev, erfcx_gradient_reference) {
  struct {
    double x;
    double d;
  } cases[] = {{-2.0, -436.89199672700738},  {-1.0, -11.14633932862008},
               {0.0, -1.1283791670955126},   {1.0, -0.27321201478389856},
               {5.0, -0.021332789764826311}, {10.0, -0.0055593122190608565}};
  for (auto c : cases) {
    stan::math::var x = c.x;
    stan::math::var y = stan::math::erfcx(x);
    y.grad();
    EXPECT_NEAR(c.d, x.adj(), 1e-12 * std::fabs(c.d));
    stan::math::recover_memory();
  }
}

// The gradient is built from the value, so it stays finite exactly where the
// value does -- this is the tail behaviour that motivates the function.
TEST_F(AgradRev, erfcx_upper_tail_gradient) {
  stan::math::var x = 30.0;
  stan::math::var y = stan::math::erfcx(x);
  y.grad();
  EXPECT_TRUE(std::isfinite(y.val()));
  EXPECT_GT(y.val(), 0.0);
  EXPECT_TRUE(std::isfinite(x.adj()));
  EXPECT_NEAR(expected_deriv(30.0), x.adj(), 1e-15);
}

TEST_F(AgradRev, erfcx_matrix) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> x(3);
  x << -1.0, 0.5, 12.0;
  stan::math::var y = stan::math::sum(stan::math::erfcx(x));
  y.grad();
  for (int i = 0; i < x.size(); ++i) {
    EXPECT_FLOAT_EQ(expected_deriv(x(i).val()), x(i).adj());
  }
}

TEST_F(AgradRev, erfcx_var_matrix) {
  Eigen::VectorXd xd(3);
  xd << -1.0, 0.5, 12.0;
  stan::math::var_value<Eigen::VectorXd> x(xd);
  auto fx = stan::math::erfcx(x);
  stan::math::var y = stan::math::sum(fx);

  double expected_sum = 0.0;
  for (int i = 0; i < xd.size(); ++i) {
    EXPECT_FLOAT_EQ(stan::math::erfcx(xd(i)), fx.val()(i));
    expected_sum += stan::math::erfcx(xd(i));
  }
  EXPECT_FLOAT_EQ(expected_sum, y.val());

  y.grad();
  for (int i = 0; i < xd.size(); ++i) {
    EXPECT_FLOAT_EQ(expected_deriv(xd(i)), x.adj()(i));
  }
}
