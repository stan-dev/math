#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <string>
#include <vector>

void test_value_and_derivatives(double expected_val, double y_dbl) {
  using stan::math::is_nan;
  using stan::math::std_normal_lccdf;
  using stan::math::var;

  std::stringstream msg_ss;
  msg_ss << "parameters: (" << y_dbl << ")";
  std::string msg = msg_ss.str();
  SCOPED_TRACE(msg);

  var y(y_dbl);
  var val = std_normal_lccdf(y);
  std::vector<var> x;
  x.push_back(y);
  std::vector<double> gradients;
  gradients.clear();
  val.grad(x, gradients);

  double e = 1e-10;
  double inv2e = 0.5 / e;
  double finite_diff
      = (std_normal_lccdf(y_dbl + e) - std_normal_lccdf(y_dbl - e)) * inv2e;

  EXPECT_FLOAT_EQ(expected_val, val.val());
  EXPECT_FALSE(is_nan(gradients[0]));
  if (!is_nan(gradients[0])) {
    if (!is_nan(finite_diff))
      EXPECT_NEAR(finite_diff, gradients[0], 1e-2);
    else
      EXPECT_FLOAT_EQ(-std::numeric_limits<double>::infinity(), gradients[0]);
  }
}

TEST(std_normal_lccdf, derivatives) {
  test_value_and_derivatives(log(0.5), 0.0);
  test_value_and_derivatives(0.0, -20.0);
  test_value_and_derivatives(0.0, -30.0);
  test_value_and_derivatives(0.0, -50.0);
  test_value_and_derivatives(-std::numeric_limits<double>::infinity(), 30.0);
}
