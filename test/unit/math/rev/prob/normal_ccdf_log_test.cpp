#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>
#include <vector>
#include <limits>
#include <string>

void test_value_and_derivatives(double expected_val, double y_dbl,
                                double mu_dbl, double sigma_dbl) {
  using stan::math::is_nan;
  using stan::math::normal_ccdf_log;
  using stan::math::var;
  std::stringstream msg_ss;
  msg_ss << "parameters: (" << y_dbl << ", " << mu_dbl << ", " << sigma_dbl
         << ")";
  std::string msg = msg_ss.str();

  SCOPED_TRACE(msg);

  var y(y_dbl);
  var mu(mu_dbl);
  var sigma(sigma_dbl);

  std::vector<double> gradients;
  var val = normal_ccdf_log(y, mu, sigma);
  std::vector<var> x;
  x.push_back(y);
  x.push_back(mu);
  x.push_back(sigma);
  gradients.clear();
  val.grad(x, gradients);

  double e = 1e-10;
  double inv2e = 0.5 / e;
  std::vector<double> finite_diffs;
  finite_diffs.resize(3);
  finite_diffs[0] = (normal_ccdf_log(y_dbl + e, mu_dbl, sigma_dbl)
                     - normal_ccdf_log(y_dbl - e, mu_dbl, sigma_dbl))
                    * inv2e;
  finite_diffs[1] = (normal_ccdf_log(y_dbl, mu_dbl + e, sigma_dbl)
                     - normal_ccdf_log(y_dbl, mu_dbl - e, sigma_dbl))
                    * inv2e;
  finite_diffs[2] = (normal_ccdf_log(y_dbl, mu_dbl, sigma_dbl + e)
                     - normal_ccdf_log(y_dbl, mu_dbl, sigma_dbl - e))
                    * inv2e;

  EXPECT_FLOAT_EQ(expected_val, val.val());
  EXPECT_FALSE(is_nan(gradients[0]));
  EXPECT_FALSE(is_nan(gradients[1]));
  EXPECT_FALSE(is_nan(gradients[2]));
  if (!is_nan(gradients[0])) {
    if (!is_nan(finite_diffs[0]))
      EXPECT_NEAR(finite_diffs[0], gradients[0], 1e-2);
    else
      EXPECT_FLOAT_EQ(-std::numeric_limits<double>::infinity(), gradients[0]);
  }
  if (!is_nan(gradients[1])) {
    if (!is_nan(finite_diffs[1]))
      EXPECT_NEAR(finite_diffs[1], gradients[1], 1e-2);
    else
      EXPECT_FLOAT_EQ(std::numeric_limits<double>::infinity(), gradients[1]);
  }
  if (!is_nan(gradients[2])) {
    if (!is_nan(finite_diffs[2]))
      EXPECT_NEAR(finite_diffs[2], gradients[2], 1e-2);
    else
      EXPECT_FLOAT_EQ(std::numeric_limits<double>::infinity(), gradients[2]);
  }
}

TEST(normal_ccdf_log, derivatives) {
  test_value_and_derivatives(log(0.5), 10.0, 10.0, 0.5);
  test_value_and_derivatives(0.0, -20.0, 10.0, 0.5);
  test_value_and_derivatives(0.0, -30.0, 10.0, 0.5);
  test_value_and_derivatives(0.0, -50.0, 10.0, 1.0);

  test_value_and_derivatives(-std::numeric_limits<double>::infinity(), 30.0,
                             10.0, 0.5);
}
