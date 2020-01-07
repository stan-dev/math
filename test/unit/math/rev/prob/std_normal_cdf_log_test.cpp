#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <string>
#include <vector>

void test_value_and_derivatives(double expected_val, double y_dbl) {
  using stan::math::is_nan;
  using stan::math::std_normal_lcdf;
  using stan::math::var;

  std::stringstream msg_ss;
  msg_ss << "for y = " << y_dbl;
  std::string msg = msg_ss.str();
  SCOPED_TRACE(msg);

  var y(y_dbl);

  std::vector<double> gradients;
  var val = std_normal_lcdf(y);
  std::vector<var> x;
  x.push_back(y);
  gradients.clear();
  val.grad(x, gradients);

  double e = 1e-10;
  double inv2e = 0.5 / e;
  double finite_diff
      = (std_normal_lcdf(y_dbl + e) - std_normal_lcdf(y_dbl - e)) * inv2e;

  EXPECT_FLOAT_EQ(expected_val, val.val());
  EXPECT_FALSE(is_nan(gradients[0]));

  // use relative error check in first case as logs of functions do not
  // necessarily converge to asymptotic values
  if (!is_nan(gradients[0])) {
    if (!is_nan(finite_diff)) {
      if (abs((finite_diff - gradients[0]) / finite_diff) > 0.01) {
        EXPECT_NEAR(finite_diff, gradients[0], 1e-2);
      }
    } else {
      EXPECT_FLOAT_EQ(std::numeric_limits<double>::infinity(), gradients[0]);
    }
  }
}

TEST(std_normal_lcdf, derivatives) {
  test_value_and_derivatives(log(0.5), 0.0);
  test_value_and_derivatives(-std::numeric_limits<double>::infinity(),
                             -1.18292e31);
  test_value_and_derivatives(0.0, 30.0);
}
