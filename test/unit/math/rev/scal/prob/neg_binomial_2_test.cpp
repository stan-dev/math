#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(ProbDistributionsNegBinomial, derivatives_not_nan) {
  using stan::math::is_nan;
  using stan::math::var;
  using stan::math::neg_binomial_2_log;

  int N = 100;
  double mu_dbl = 8;
  double phi_dbl = 1;
  
  var mu(mu_dbl);
  var phi(phi_dbl);

  for (int k = 0; k < 20; ++k) {
    var val = neg_binomial_2_log(N, mu, phi);

    std::vector<var> x;
    x.push_back(mu);
    x.push_back(phi);
    
    std::vector<double> gradients;
    val.grad(x, gradients);

    EXPECT_FALSE(is_nan(gradients[0]));
    EXPECT_FALSE(is_nan(gradients[1]));

    phi *= 10;
  }

}

void test_value_and_derivatives(double expected_val,
                                double y_dbl, double mu_dbl, double sigma_dbl) {
  using stan::math::is_nan;
  using stan::math::var;  
  using stan::math::normal_cdf;
  std::stringstream msg_ss;
  msg_ss << "parameters: (" << y_dbl << ", " << mu_dbl << ", " << sigma_dbl << ")";
  std::string msg = msg_ss.str();

  SCOPED_TRACE(msg);
  
  var y(y_dbl);
  var mu(mu_dbl);
  var sigma(sigma_dbl);

  std::vector<double> gradients;
  var val = normal_cdf(y, mu, sigma);
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
  finite_diffs[0] = (normal_cdf(y_dbl + e, mu_dbl, sigma_dbl)
                     - normal_cdf(y_dbl - e, mu_dbl, sigma_dbl)) * inv2e;
  finite_diffs[1] = (normal_cdf(y_dbl, mu_dbl + e, sigma_dbl)
                     - normal_cdf(y_dbl, mu_dbl - e, sigma_dbl)) * inv2e;
  finite_diffs[2] = (normal_cdf(y_dbl, mu_dbl, sigma_dbl + e)
                     - normal_cdf(y_dbl, mu_dbl, sigma_dbl - e)) * inv2e;

  
  EXPECT_FLOAT_EQ(expected_val, val.val());
  EXPECT_FALSE(is_nan(gradients[0]));
  EXPECT_FALSE(is_nan(gradients[1]));
  EXPECT_FALSE(is_nan(gradients[2]));
  if (!is_nan(gradients[0])) {
    if (!is_nan(finite_diffs[0]))
      EXPECT_NEAR(finite_diffs[0], gradients[0], 1e-2);
    else
      EXPECT_FLOAT_EQ(0.0, gradients[0]);
  }
  if (!is_nan(gradients[1])) {
    if (!is_nan(finite_diffs[1]))
      EXPECT_NEAR(finite_diffs[1], gradients[1], 1e-2);
    else
      EXPECT_FLOAT_EQ(0.0, gradients[1]);
  }
  if (!is_nan(gradients[2])) {
    if (!is_nan(finite_diffs[2]))
      EXPECT_NEAR(finite_diffs[2], gradients[2], 1e-2);
    else
      EXPECT_FLOAT_EQ(0.0, gradients[2]);
  }
}


// TEST(normal_cdf, derivatives) {
//   test_value_and_derivatives(0.5, 10.0, 10.0, 0.5);
//   test_value_and_derivatives(0.0, -20.0, 10.0, 0.5);
//   test_value_and_derivatives(0.0, -30.0, 10.0, 0.5);
//   test_value_and_derivatives(0.0, -50.0, 10.0, 1.0);
// 
//   test_value_and_derivatives(1.0, 30.0, 10.0, 0.5);
// }
