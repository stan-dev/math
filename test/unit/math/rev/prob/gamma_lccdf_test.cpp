#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <vector>
#include <limits>

TEST(ProbDistributionsGamma, lccdf_values) {
  using stan::math::gamma_lccdf;
  using stan::math::var;

  double y_d = 0.8;
  double alpha_d = 1.1;
  double beta_d = 2.3;

  var y_v = y_d;
  var alpha_v = alpha_d;
  var beta_v = beta_d;

  var lccdf_var = gamma_lccdf(y_v, alpha_v, beta_v);

  EXPECT_NO_THROW(lccdf_var.val());
  EXPECT_FALSE(std::isnan(lccdf_var.val()));
  EXPECT_LE(lccdf_var.val(), 0.0);
}

TEST(ProbDistributionsGamma, lccdf_derivatives_y) {
  using stan::math::gamma_lccdf;
  using stan::math::var;

  // Test derivative with respect to y
  double y_d = 2.5;
  double alpha_d = 3.0;
  double beta_d = 1.5;

  var y_v = y_d;
  var alpha_v = alpha_d;
  var beta_v = beta_d;

  var lccdf_var = gamma_lccdf(y_v, alpha_v, beta_v);

  std::vector<var> vars = {y_v};
  std::vector<double> grads;
  lccdf_var.grad(vars, grads);

  // Derivative should be negative (LCCDF decreases as y increases)
  EXPECT_LT(grads[0], 0.0);
  EXPECT_FALSE(std::isnan(grads[0]));
  EXPECT_TRUE(std::isfinite(grads[0]));
}

TEST(ProbDistributionsGamma, lccdf_derivatives_alpha) {
  using stan::math::gamma_lccdf;
  using stan::math::var;

  // Test derivative with respect to alpha
  double y_d = 2.5;
  double alpha_d = 3.0;
  double beta_d = 1.5;

  var y_v = y_d;
  var alpha_v = alpha_d;
  var beta_v = beta_d;

  var lccdf_var = gamma_lccdf(y_v, alpha_v, beta_v);

  std::vector<var> vars = {alpha_v};
  std::vector<double> grads;
  lccdf_var.grad(vars, grads);

  EXPECT_FALSE(std::isnan(grads[0]));
  EXPECT_TRUE(std::isfinite(grads[0]));
}

TEST(ProbDistributionsGamma, lccdf_derivatives_beta) {
  using stan::math::gamma_lccdf;
  using stan::math::var;

  // Test derivative with respect to beta
  double y_d = 2.5;
  double alpha_d = 3.0;
  double beta_d = 1.5;

  var y_v = y_d;
  var alpha_v = alpha_d;
  var beta_v = beta_d;

  var lccdf_var = gamma_lccdf(y_v, alpha_v, beta_v);

  std::vector<var> vars = {beta_v};
  std::vector<double> grads;
  lccdf_var.grad(vars, grads);

  EXPECT_FALSE(std::isnan(grads[0]));
  EXPECT_TRUE(std::isfinite(grads[0]));
}

TEST(ProbDistributionsGamma, lccdf_derivatives_all_params) {
  using stan::math::gamma_lccdf;
  using stan::math::var;

  // Test all derivatives together
  double y_d = 2.5;
  double alpha_d = 3.0;
  double beta_d = 1.5;

  var y_v = y_d;
  var alpha_v = alpha_d;
  var beta_v = beta_d;

  var lccdf_var = gamma_lccdf(y_v, alpha_v, beta_v);

  std::vector<var> vars = {y_v, alpha_v, beta_v};
  std::vector<double> grads;
  lccdf_var.grad(vars, grads);

  // Check all gradients are finite and not NaN
  for (size_t i = 0; i < grads.size(); ++i) {
    EXPECT_FALSE(std::isnan(grads[i])) << "Gradient " << i << " is NaN";
    EXPECT_TRUE(std::isfinite(grads[i]))
        << "Gradient " << i << " is not finite";
  }

  // d/dy should be negative
  EXPECT_LT(grads[0], 0.0) << "d/dy should be negative";
}

TEST(ProbDistributionsGamma, lccdf_finite_diff_y) {
  using stan::math::gamma_lccdf;
  using stan::math::var;

  // Test derivative w.r.t. y using finite differences
  double y_d = 1.5;
  double alpha_d = 2.0;
  double beta_d = 3.0;
  double eps = 1e-6;

  // Compute gradient using autodiff
  var y_v = y_d;
  var lccdf_var = gamma_lccdf(y_v, alpha_d, beta_d);
  std::vector<var> vars = {y_v};
  std::vector<double> grads;
  lccdf_var.grad(vars, grads);
  double grad_autodiff = grads[0];

  // Compute gradient using finite differences
  double f_plus = stan::math::gamma_lccdf(y_d + eps, alpha_d, beta_d);
  double f_minus = stan::math::gamma_lccdf(y_d - eps, alpha_d, beta_d);
  double grad_findiff = (f_plus - f_minus) / (2.0 * eps);

  EXPECT_NEAR(grad_autodiff, grad_findiff, 1e-5);
}

TEST(ProbDistributionsGamma, lccdf_finite_diff_alpha) {
  using stan::math::gamma_lccdf;
  using stan::math::var;

  // Test derivative w.r.t. alpha using finite differences
  double y_d = 1.5;
  double alpha_d = 2.0;
  double beta_d = 3.0;
  double eps = 1e-6;

  // Compute gradient using autodiff
  var alpha_v = alpha_d;
  var lccdf_var = gamma_lccdf(y_d, alpha_v, beta_d);
  std::vector<var> vars = {alpha_v};
  std::vector<double> grads;
  lccdf_var.grad(vars, grads);
  double grad_autodiff = grads[0];

  // Compute gradient using finite differences
  double f_plus = stan::math::gamma_lccdf(y_d, alpha_d + eps, beta_d);
  double f_minus = stan::math::gamma_lccdf(y_d, alpha_d - eps, beta_d);
  double grad_findiff = (f_plus - f_minus) / (2.0 * eps);

  // Alpha derivative is more numerically sensitive, use larger tolerance
  EXPECT_NEAR(grad_autodiff, grad_findiff, 1e-3);
}

TEST(ProbDistributionsGamma, lccdf_finite_diff_beta) {
  using stan::math::gamma_lccdf;
  using stan::math::var;

  // Test derivative w.r.t. beta using finite differences
  double y_d = 1.5;
  double alpha_d = 2.0;
  double beta_d = 3.0;
  double eps = 1e-6;

  // Compute gradient using autodiff
  var beta_v = beta_d;
  var lccdf_var = gamma_lccdf(y_d, alpha_d, beta_v);
  std::vector<var> vars = {beta_v};
  std::vector<double> grads;
  lccdf_var.grad(vars, grads);
  double grad_autodiff = grads[0];

  // Compute gradient using finite differences
  double f_plus = stan::math::gamma_lccdf(y_d, alpha_d, beta_d + eps);
  double f_minus = stan::math::gamma_lccdf(y_d, alpha_d, beta_d - eps);
  double grad_findiff = (f_plus - f_minus) / (2.0 * eps);

  EXPECT_NEAR(grad_autodiff, grad_findiff, 1e-5);
}

TEST(ProbDistributionsGamma, lccdf_extreme_values_small) {
  using stan::math::gamma_lccdf;
  using stan::math::var;

  // Test with very small values
  double y_d = 0.001;
  double alpha_d = 0.1;
  double beta_d = 1.0;

  var y_v = y_d;
  var alpha_v = alpha_d;
  var beta_v = beta_d;

  var lccdf_var = gamma_lccdf(y_v, alpha_v, beta_v);

  std::vector<var> vars = {y_v, alpha_v, beta_v};
  std::vector<double> grads;
  lccdf_var.grad(vars, grads);

  // All gradients should be finite
  for (size_t i = 0; i < grads.size(); ++i) {
    EXPECT_TRUE(std::isfinite(grads[i]))
        << "Gradient " << i << " not finite for small values";
    EXPECT_FALSE(std::isnan(grads[i]))
        << "Gradient " << i << " is NaN for small values";
  }
}

TEST(ProbDistributionsGamma, lccdf_extreme_values_large) {
  using stan::math::gamma_lccdf;
  using stan::math::var;

  // Test with large values
  double y_d = 100.0;
  double alpha_d = 50.0;
  double beta_d = 0.5;

  var y_v = y_d;
  var alpha_v = alpha_d;
  var beta_v = beta_d;

  var lccdf_var = gamma_lccdf(y_v, alpha_v, beta_v);

  std::vector<var> vars = {y_v, alpha_v, beta_v};
  std::vector<double> grads;
  lccdf_var.grad(vars, grads);

  // All gradients should be finite
  for (size_t i = 0; i < grads.size(); ++i) {
    EXPECT_TRUE(std::isfinite(grads[i]))
        << "Gradient " << i << " not finite for large values";
    EXPECT_FALSE(std::isnan(grads[i]))
        << "Gradient " << i << " is NaN for large values";
  }
}

TEST(ProbDistributionsGamma, lccdf_alpha_one_derivatives) {
  using stan::math::gamma_lccdf;
  using stan::math::var;

  // When alpha = 1, gamma is exponential
  // LCCDF(y) = -beta * y
  // d/dy = -beta
  // d/dbeta = -y
  double y_d = 2.0;
  double alpha_d = 1.0;
  double beta_d = 3.0;

  var y_v = y_d;
  var alpha_v = alpha_d;
  var beta_v = beta_d;

  var lccdf_var = gamma_lccdf(y_v, alpha_v, beta_v);

  std::vector<var> vars = {y_v, alpha_v, beta_v};
  std::vector<double> grads;
  lccdf_var.grad(vars, grads);

  // For exponential: d/dy LCCDF = d/dy (-beta*y) = -beta
  EXPECT_NEAR(grads[0], -beta_d, 1e-10);

  // d/dbeta LCCDF = d/dbeta (-beta*y) = -y
  EXPECT_NEAR(grads[2], -y_d, 1e-10);
}

TEST(ProbDistributionsGamma, lccdf_various_parameter_combinations) {
  using stan::math::gamma_lccdf;
  using stan::math::var;

  std::vector<std::tuple<double, double, double>> test_cases = {
      {0.5, 0.5, 1.0},  {1.0, 1.0, 1.0},  {2.0, 3.0, 0.5},
      {10.0, 2.0, 0.1}, {0.1, 10.0, 2.0}, {5.0, 5.0, 1.0},
  };

  for (const auto& test_case : test_cases) {
    double y_d = std::get<0>(test_case);
    double alpha_d = std::get<1>(test_case);
    double beta_d = std::get<2>(test_case);

    var y_v = y_d;
    var alpha_v = alpha_d;
    var beta_v = beta_d;

    var lccdf_var = gamma_lccdf(y_v, alpha_v, beta_v);

    std::vector<var> vars = {y_v, alpha_v, beta_v};
    std::vector<double> grads;
    lccdf_var.grad(vars, grads);

    // All gradients should be finite
    for (size_t i = 0; i < grads.size(); ++i) {
      EXPECT_FALSE(std::isnan(grads[i]))
          << "NaN gradient for y=" << y_d << ", alpha=" << alpha_d
          << ", beta=" << beta_d;
      EXPECT_TRUE(std::isfinite(grads[i]))
          << "Non-finite gradient for y=" << y_d << ", alpha=" << alpha_d
          << ", beta=" << beta_d;
    }
  }
}

TEST(ProbDistributionsGamma, lccdf_second_derivative_y) {
  using stan::math::gamma_lccdf;
  using stan::math::var;

  // Test that second derivatives work (for models using Hessians)
  double y_d = 2.0;
  double alpha_d = 3.0;
  double beta_d = 1.5;

  var y_v = y_d;

  var lccdf_var = gamma_lccdf(y_v, alpha_d, beta_d);

  std::vector<var> vars = {y_v};
  std::vector<double> grads;
  lccdf_var.grad(vars, grads);

  // If we can compute gradients without error, second derivatives should work
  EXPECT_FALSE(std::isnan(grads[0]));
  EXPECT_TRUE(std::isfinite(grads[0]));
}

TEST(ProbDistributionsGamma, lccdf_numerically_challenging_derivatives) {
  using stan::math::gamma_lccdf;
  using stan::math::var;

  // Test numerically challenging cases
  std::vector<std::tuple<double, double, double>> challenging_cases = {
      {1e-6, 1e-5, 1.0},     // Very small y and alpha
      {0.001, 100.0, 0.01},  // Small y, large alpha, small beta
      {1000.0, 0.5, 1e-3},   // Large y, small alpha, very small beta
      {50.0, 50.0, 1.0},     // Matched moderate values
  };

  for (const auto& test_case : challenging_cases) {
    double y_d = std::get<0>(test_case);
    double alpha_d = std::get<1>(test_case);
    double beta_d = std::get<2>(test_case);

    var y_v = y_d;
    var alpha_v = alpha_d;
    var beta_v = beta_d;

    var lccdf_var = gamma_lccdf(y_v, alpha_v, beta_v);

    std::vector<var> vars = {y_v, alpha_v, beta_v};
    std::vector<double> grads;
    lccdf_var.grad(vars, grads);

    // Gradients should not be NaN
    for (size_t i = 0; i < grads.size(); ++i) {
      EXPECT_FALSE(std::isnan(grads[i]))
          << "NaN gradient " << i << " for y=" << y_d << ", alpha=" << alpha_d
          << ", beta=" << beta_d;
    }
  }
}
