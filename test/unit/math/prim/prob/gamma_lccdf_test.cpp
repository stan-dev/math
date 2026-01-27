#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(ProbGamma, lccdf_works) {
  using stan::math::gamma_lccdf;

  double y = 0.8;
  double alpha = 1.1;
  double beta = 2.3;

  EXPECT_NO_THROW(gamma_lccdf(y, alpha, beta));
}

TEST(ProbGamma, lccdf_zero_y) {
  using stan::math::gamma_lccdf;

  // When y = 0, LCCDF(0) = log(P(Y > 0)) = log(1) = 0
  // For continuous distribution, P(Y > 0) = 1
  double alpha = 1.5;
  double beta = 2.0;

  double result = gamma_lccdf(0.0, alpha, beta);
  EXPECT_EQ(result, 0.0);
}

TEST(ProbGamma, lccdf_large_y) {
  using stan::math::gamma_lccdf;

  // When y is very large, CDF approaches 1, so LCCDF = log(1-1) = log(0) = -inf
  double alpha = 1.5;
  double beta = 2.0;
  double y = 1e10;

  double result = gamma_lccdf(y, alpha, beta);

  // Should be a very large negative number (approaching -infinity)
  EXPECT_LT(result, -1000.0);
}

TEST(ProbGamma, lccdf_infinity_y) {
  using stan::math::gamma_lccdf;
  using stan::math::negative_infinity;

  // When y = infinity, LCCDF = log(P(Y > ∞)) = log(0) = -∞
  double alpha = 1.5;
  double beta = 2.0;
  double y = std::numeric_limits<double>::infinity();

  double result = gamma_lccdf(y, alpha, beta);
  EXPECT_EQ(result, negative_infinity());
}

TEST(ProbGamma, lccdf_small_alpha_small_y) {
  using stan::math::gamma_lccdf;

  // Small alpha, small y - numerically challenging
  double y = 0.001;
  double alpha = 0.1;
  double beta = 1.0;

  double result = gamma_lccdf(y, alpha, beta);

  // Should be finite and negative
  EXPECT_TRUE(std::isfinite(result));
  EXPECT_LT(result, 0.0);
}

TEST(ProbGamma, lccdf_large_alpha_large_y) {
  using stan::math::gamma_lccdf;

  // Large alpha, large y
  double y = 100.0;
  double alpha = 50.0;
  double beta = 0.5;

  double result = gamma_lccdf(y, alpha, beta);

  // Should be finite
  EXPECT_TRUE(std::isfinite(result));
}

TEST(ProbGamma, lccdf_alpha_one) {
  using stan::math::gamma_lccdf;
  using std::exp;
  using std::log;

  // When alpha = 1, gamma becomes exponential
  // For exponential with rate beta: LCCDF(y) = log(1 - (1-exp(-beta*y))) =
  // log(exp(-beta*y)) = -beta*y
  double y = 2.0;
  double alpha = 1.0;
  double beta = 3.0;

  double result = gamma_lccdf(y, alpha, beta);
  double expected = -beta * y;  // = -6.0

  EXPECT_NEAR(result, expected, 1e-10);
}

TEST(ProbGamma, lccdf_various_values) {
  using stan::math::gamma_lccdf;

  // Test a variety of parameter combinations
  std::vector<std::tuple<double, double, double>> test_cases = {
      {0.5, 0.5, 1.0},     // Small y, small alpha
      {1.0, 1.0, 1.0},     // All ones
      {2.0, 3.0, 0.5},     // Moderate values
      {10.0, 2.0, 0.1},    // Large y, small beta
      {0.1, 10.0, 2.0},    // Small y, large alpha
      {5.0, 5.0, 1.0},     // Equal alpha and y
      {0.01, 0.5, 10.0},   // Small y, large beta
      {100.0, 100.0, 1.0}  // Large matched values
  };

  for (const auto& test_case : test_cases) {
    double y = std::get<0>(test_case);
    double alpha = std::get<1>(test_case);
    double beta = std::get<2>(test_case);

    double result = gamma_lccdf(y, alpha, beta);

    // All results should be finite and <= 0
    EXPECT_TRUE(std::isfinite(result))
        << "Failed for y=" << y << ", alpha=" << alpha << ", beta=" << beta;
    EXPECT_LE(result, 0.0) << "Failed for y=" << y << ", alpha=" << alpha
                           << ", beta=" << beta;
  }
}

TEST(ProbGamma, lccdf_extreme_small_values) {
  using stan::math::gamma_lccdf;

  // Very small but non-zero values
  double y = 1e-10;
  double alpha = 1e-5;
  double beta = 1.0;

  double result = gamma_lccdf(y, alpha, beta);

  EXPECT_TRUE(std::isfinite(result));
}

TEST(ProbGamma, lccdf_extreme_large_alpha) {
  using stan::math::gamma_lccdf;

  // Very large alpha (approaches normal distribution)
  double y = 1000.0;
  double alpha = 1000.0;
  double beta = 1.0;

  double result = gamma_lccdf(y, alpha, beta);

  EXPECT_TRUE(std::isfinite(result));
}

TEST(ProbGamma, lccdf_monotonic_in_y) {
  using stan::math::gamma_lccdf;

  // LCCDF should be monotonically decreasing in y
  double alpha = 2.0;
  double beta = 1.5;

  double y1 = 1.0;
  double y2 = 2.0;
  double y3 = 3.0;

  double lccdf1 = gamma_lccdf(y1, alpha, beta);
  double lccdf2 = gamma_lccdf(y2, alpha, beta);
  double lccdf3 = gamma_lccdf(y3, alpha, beta);

  EXPECT_GT(lccdf1, lccdf2);
  EXPECT_GT(lccdf2, lccdf3);
}

TEST(ProbGamma, lccdf_consistency_with_cdf) {
  using stan::math::gamma_cdf;
  using stan::math::gamma_lccdf;
  using std::log;

  // Test that lccdf(y) ≈ log(1 - cdf(y))
  double y = 1.5;
  double alpha = 2.5;
  double beta = 1.8;

  double lccdf_val = gamma_lccdf(y, alpha, beta);
  double cdf_val = gamma_cdf(y, alpha, beta);
  double expected = log(1.0 - cdf_val);

  EXPECT_NEAR(lccdf_val, expected, 1e-10);
}

TEST(ProbGamma, lccdf_numerically_challenging) {
  using stan::math::gamma_lccdf;

  // Test cases that might cause numerical issues
  std::vector<std::tuple<double, double, double>> challenging_cases = {
      {1e-8, 1e-6, 1.0},      // Very small y and alpha
      {1e-6, 100.0, 1e-3},    // Very small y, large alpha, small beta
      {1000.0, 0.1, 1e-4},    // Large y, small alpha, very small beta
      {50.0, 50.0, 1.0},      // Matched moderate values
      {0.001, 0.001, 100.0},  // Small y and alpha, large beta
      {1e6, 10.0, 1e-6},      // Very large y, moderate alpha, very small beta
  };

  for (const auto& test_case : challenging_cases) {
    double y = std::get<0>(test_case);
    double alpha = std::get<1>(test_case);
    double beta = std::get<2>(test_case);

    double result = gamma_lccdf(y, alpha, beta);

    // Should not be NaN
    EXPECT_FALSE(std::isnan(result))
        << "NaN for y=" << y << ", alpha=" << alpha << ", beta=" << beta;

    // Should be <= 0 (log of probability)
    EXPECT_LE(result, 0.0) << "Positive value for y=" << y
                           << ", alpha=" << alpha << ", beta=" << beta;
  }
}

TEST(ProbGamma, lccdf_shape_zero_throws) {
  using stan::math::gamma_lccdf;

  // alpha (shape) must be positive
  EXPECT_THROW(gamma_lccdf(1.0, 0.0, 1.0), std::domain_error);
  EXPECT_THROW(gamma_lccdf(1.0, -1.0, 1.0), std::domain_error);
}

TEST(ProbGamma, lccdf_rate_zero_throws) {
  using stan::math::gamma_lccdf;

  // beta (rate) must be positive
  EXPECT_THROW(gamma_lccdf(1.0, 1.0, 0.0), std::domain_error);
  EXPECT_THROW(gamma_lccdf(1.0, 1.0, -1.0), std::domain_error);
}

TEST(ProbGamma, lccdf_negative_y_throws) {
  using stan::math::gamma_lccdf;

  // y must be non-negative
  EXPECT_THROW(gamma_lccdf(-1.0, 1.0, 1.0), std::domain_error);
  EXPECT_THROW(gamma_lccdf(-0.001, 1.0, 1.0), std::domain_error);
}
