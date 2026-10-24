#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <boost/math/interpolators/cubic_hermite.hpp>
#include <boost/math/interpolators/makima.hpp>
#include <limits>
#include <tuple>
#include <cmath>

namespace {
// Standard test data with diverse characteristics
const Eigen::VectorXd kTestX = (Eigen::VectorXd(4) << 0.0, 1.0, 2.5, 4.0).finished();
const Eigen::VectorXd kTestY = (Eigen::VectorXd(4) << 10.0, 12.0, 9.0, 11.0).finished();
const Eigen::VectorXd kTestDydx = (Eigen::VectorXd(4) << 1.0, 0.5, -0.5, 1.5).finished();

// Test data with very close spacing for numerical stability (1e-6 intervals)
const Eigen::VectorXd kCloseSpaceX = (Eigen::VectorXd(3) << 0.0, 1e-6, 1.0).finished();
const Eigen::VectorXd kCloseSpaceY = (Eigen::VectorXd(3) << 5.0, 5.000001, 10.0).finished();
const Eigen::VectorXd kCloseSpaceDydx = (Eigen::VectorXd(3) << 2.0, 1.5, 3.0).finished();

/**
 * Helper function to create a Boost cubic Hermite interpolator from Eigen vectors.
 * 
 * The Boost interpolator serves as the ground truth reference implementation.
 * It uses the same cubic Hermite basis functions as the Stan implementation.
 */
auto create_boost_cubic_hermite_interpolator(const Eigen::VectorXd& x, const Eigen::VectorXd& y, 
                               const Eigen::VectorXd& dydx) {
  std::vector<double> x_vec(x.data(), x.data() + x.size());
  std::vector<double> y_vec(y.data(), y.data() + y.size());
  std::vector<double> dydx_vec(dydx.data(), dydx.data() + dydx.size());
  
  return boost::math::interpolators::cubic_hermite<std::vector<double>>(
    std::move(x_vec), std::move(y_vec), std::move(dydx_vec));
}
}  // namespace

// ============================================================================
// SECTION 1: Setup Function Tests
// ============================================================================

/**
 * TEST: interpolate1DCubicHermiteSetupStorage
 * 
 * Verifies that interpolate_1d_cubic_hermite_setup correctly stores and validates
 * the input knot locations, values, and derivatives in the output tuple.
 * 
 * TUPLE STRUCTURE: (xk, yk, dydxk)
 * 
 * The setup function should:
 * 1. Validate all inputs (size matching, finiteness, ordering)
 * 2. Store xk, yk, and dydxk in the output tuple
 * 3. Preserve the original values
 */
TEST(MathPrimFun, interpolate1DCubicHermiteSetupStorage) {
  auto W = stan::math::interpolate_1d_cubic_hermite_setup(kTestX, kTestY, kTestDydx);
  
  // Extract components from the tuple structure
  const auto& xk = std::get<0>(W);  // Knot locations
  const auto& yk = std::get<1>(W);  // Knot values
  const auto& dydxk = std::get<2>(W);  // Derivative values

  // Verify sizes
  EXPECT_EQ(kTestX.size(), xk.size());
  EXPECT_EQ(kTestY.size(), yk.size());
  EXPECT_EQ(kTestDydx.size(), dydxk.size());

  // Verify knot locations are preserved
  for (Eigen::Index i = 0; i < kTestX.size(); ++i) {
    EXPECT_FLOAT_EQ(kTestX(i), xk(i));
  }

  // Verify knot values are preserved
  for (Eigen::Index i = 0; i < kTestY.size(); ++i) {
    EXPECT_FLOAT_EQ(kTestY(i), yk(i));
  }

  // Verify derivative values are preserved
  for (Eigen::Index i = 0; i < kTestDydx.size(); ++i) {
    EXPECT_FLOAT_EQ(kTestDydx(i), dydxk(i));
  }
}



// ============================================================================
// SECTION 2: Evaluation Function Tests (using Boost as Ground Truth)
// ============================================================================

/**
 * TEST: interpolate1DCubicHermiteEval
 * 
 * Verifies that interpolate_1d_cubic_spline correctly evaluates the cubic
 * Hermite interpolant at various locations using the standard Hermite basis
 * functions:
 *   y(x) = (1-t)^2 * (y0*(1+2*t) + h*s0*t) + t^2 * (y1*(3-2*t) + h*s1*(t-1))
 * where t = (x - xk(i)) / h.
 * 
 * IMPORTANT: Ground truth values are computed using the Boost cubic Hermite
 * interpolator as a reference implementation. Both Stan and Boost use the same
 * cubic Hermite basis functions, so they should produce identical results.
 * 
 * Test data:
 *   xk = [0.0, 1.0, 2.5, 4.0]
 *   yk = [10.0, 12.0, 9.0, 11.0]
 *   dydxk = [1.0, 0.5, -0.5, 1.5]
 */
TEST(MathPrimFun, interpolate1DCubicHermiteEval) {
  auto W = stan::math::interpolate_1d_cubic_hermite_setup(kTestX, kTestY, kTestDydx);
  auto boost_interp = create_boost_cubic_hermite_interpolator(kTestX, kTestY, kTestDydx);

  // Test 1: Endpoint evaluation at x = xk(0) = 0.0
  {
    double x = 0.0;
    double expected = boost_interp(x);
    double computed = stan::math::interpolate_1d_cubic_spline(x, W);
    EXPECT_NEAR(expected, computed, 1e-10);
  }

  // Test 2: Endpoint evaluation at x = xk(3) = 4.0
  {
    double x = 4.0;
    double expected = boost_interp(x);
    double computed = stan::math::interpolate_1d_cubic_spline(x, W);
    EXPECT_NEAR(expected, computed, 1e-10);
  }

  // Test 3: Intermediate point in interval [0.0, 1.0]
  {
    double x = 0.5;
    double expected = boost_interp(x);
    double computed = stan::math::interpolate_1d_cubic_spline(x, W);
    EXPECT_NEAR(expected, computed, 1e-10);
  }

  // Test 4: Intermediate point in interval [1.0, 2.5]
  {
    double x = 1.5;
    double expected = boost_interp(x);
    double computed = stan::math::interpolate_1d_cubic_spline(x, W);
    EXPECT_NEAR(expected, computed, 1e-10);
  }

  // Test 5: Intermediate point in interval [2.5, 4.0]
  {
    double x = 3.0;
    double expected = boost_interp(x);
    double computed = stan::math::interpolate_1d_cubic_spline(x, W);
    EXPECT_NEAR(expected, computed, 1e-10);
  }

  // Test 6: Point very close to first knot
  {
    double x = 1e-6;
    double expected = boost_interp(x);
    double computed = stan::math::interpolate_1d_cubic_spline(x, W);
    EXPECT_NEAR(expected, computed, 1e-10);
  }

  // Test 7: Point very close to last knot
  {
    double x = 4-1e-6;
    double expected = boost_interp(x);
    double computed = stan::math::interpolate_1d_cubic_spline(x, W);
    EXPECT_NEAR(expected, computed, 1e-10);
  }
}



/**
 * TEST: interpolate1DCubicHermiteEvalCloseSpacing
 * 
 * Tests evaluation with 1e-6 interval spacing for numerical stability.
 * Verifies robustness at precision limits.
 */
TEST(MathPrimFun, interpolate1DCubicHermiteEvalCloseSpacing) {
  auto W = stan::math::interpolate_1d_cubic_hermite_setup(kCloseSpaceX, kCloseSpaceY, kCloseSpaceDydx);
  auto boost_interp = create_boost_cubic_hermite_interpolator(kCloseSpaceX, kCloseSpaceY, kCloseSpaceDydx);

  // Test at knots and intermediate point
  double test_points[] = {0.0, 5e-7, 1e-6, 1.0};
  for (double x : test_points) {
    double expected = boost_interp(x);
    double computed = stan::math::interpolate_1d_cubic_spline(x, W);
    EXPECT_NEAR(expected, computed, 1e-8) << "Failed at x = " << x;
  }
}

/**
 * TEST: interpolate1DCubicHermiteEvalMinimalKnots
 * 
 * Tests evaluation with exactly 2 knots (the minimum).
 */
TEST(MathPrimFun, interpolate1DCubicHermiteEvalMinimalKnots) {
  Eigen::VectorXd x = (Eigen::VectorXd(2) << 0.0, 1.0).finished();
  Eigen::VectorXd y = (Eigen::VectorXd(2) << 10.0, 12.0).finished();
  Eigen::VectorXd dydx = (Eigen::VectorXd(2) << 1.0, 0.5).finished();

  auto W = stan::math::interpolate_1d_cubic_hermite_setup(x, y, dydx);
  auto boost_interp = create_boost_cubic_hermite_interpolator(x, y, dydx);

  // Test at endpoints and midpoint
  double test_points[] = {0.0, 0.5, 1.0};
  for (double x_val : test_points) {
    double expected = boost_interp(x_val);
    double computed = stan::math::interpolate_1d_cubic_spline(x_val, W);
    EXPECT_NEAR(expected, computed, 1e-10);
  }
}

/**
 * TEST: interpolate1DCubicHermiteEvalNegativeDomainAndDerivatives
 * 
 * Tests evaluation with negative x values and negative derivatives.
 */
TEST(MathPrimFun, interpolate1DCubicHermiteEvalNegativeDomainAndDerivatives) {
  Eigen::VectorXd x = (Eigen::VectorXd(3) << -2.0, 0.0, 2.0).finished();
  Eigen::VectorXd y = (Eigen::VectorXd(3) << 5.0, 10.0, 8.0).finished();
  Eigen::VectorXd dydx = (Eigen::VectorXd(3) << -1.5, 0.0, -0.5).finished();

  auto W = stan::math::interpolate_1d_cubic_hermite_setup(x, y, dydx);
  auto boost_interp = create_boost_cubic_hermite_interpolator(x, y, dydx);

  std::vector<double> test_points = {-2.0, -1.0, 0.0, 1.0, 2.0};
  
  for (double x_val : test_points) {
    double expected = boost_interp(x_val);
    double computed = stan::math::interpolate_1d_cubic_spline(x_val, W);
    EXPECT_NEAR(expected, computed, 1e-10);
  }
}

/**
 * TEST: interpolate1DCubicHermiteEvalLargeDerivatives
 * 
 * Tests evaluation with large derivative values to verify numerical stability.
 */
TEST(MathPrimFun, interpolate1DCubicHermiteEvalLargeDerivatives) {
  Eigen::VectorXd x = (Eigen::VectorXd(3) << 0.0, 1.0, 2.0).finished();
  Eigen::VectorXd y = (Eigen::VectorXd(3) << 10.0, 15.0, 12.0).finished();
  Eigen::VectorXd dydx = (Eigen::VectorXd(3) << 100.0, -50.0, 75.0).finished();

  auto W = stan::math::interpolate_1d_cubic_hermite_setup(x, y, dydx);
  auto boost_interp = create_boost_cubic_hermite_interpolator(x, y, dydx);

  double x_val = 0.5;
  double expected = boost_interp(x_val);
  double computed = stan::math::interpolate_1d_cubic_spline(x_val, W);
  EXPECT_NEAR(expected, computed, 1e-9);
}

// ============================================================================
// SECTION 3: Knot Value and Derivative Boundary Tests
// ============================================================================

/**
 * TEST: interpolate1DCubicHermiteEvalAtKnots
 * 
 * Verifies that the interpolant passes through the knot values.
 * This is a fundamental property of Hermite interpolation: the interpolant
 * must pass through all specified knot points.
 */
TEST(MathPrimFun, interpolate1DCubicHermiteEvalAtKnots) {
  auto W = stan::math::interpolate_1d_cubic_hermite_setup(kTestX, kTestY, kTestDydx);

  // Verify that evaluation at each knot returns the knot value
  for (Eigen::Index i = 0; i < kTestX.size(); ++i) {
    double expected_y = kTestY(i);
    double computed_y = stan::math::interpolate_1d_cubic_spline(kTestX(i), W);
    EXPECT_NEAR(expected_y, computed_y, 1e-10) << "Failed at knot " << i;
  }
}



// ============================================================================
// SECTION 4: Exception Handling Tests for Setup Function
// ============================================================================

/**
 * TEST: interpolate1DCubicHermiteSetupThrows
 * 
 * Verifies that interpolate_1d_cubic_hermite_setup throws appropriate exceptions
 * for invalid inputs.
 */
TEST(MathPrimFun, interpolate1DCubicHermiteSetupThrowsInvalidArgument) {
  // Size mismatch: xk and yk have different sizes
  {
    Eigen::VectorXd x_bad = (Eigen::VectorXd(3) << 0.0, 1.0, 2.5).finished();
    Eigen::VectorXd y_bad = (Eigen::VectorXd(4) << 10.0, 12.0, 9.0, 11.0).finished();
    Eigen::VectorXd dydx_good = kTestDydx;

    EXPECT_THROW(stan::math::interpolate_1d_cubic_hermite_setup(x_bad, y_bad, dydx_good),
                 std::invalid_argument);
  }

  // Size mismatch: xk and dydxk have different sizes
  {
    Eigen::VectorXd x_good = kTestX;
    Eigen::VectorXd y_good = kTestY;
    Eigen::VectorXd dydx_bad = (Eigen::VectorXd(3) << 1.0, 0.5, -0.5).finished();

    EXPECT_THROW(stan::math::interpolate_1d_cubic_hermite_setup(x_good, y_good, dydx_bad),
                 std::invalid_argument);
  }
}

/**
 * TEST: interpolate1DCubicHermiteSetupThrowsDomainErrors
 * 
 * Verifies that interpolate_1d_cubic_hermite_setup throws std::domain_error
 * for invalid domain inputs.
 */
TEST(MathPrimFun, interpolate1DCubicHermiteSetupThrowsDomainErrors) {
  // Not enough knots
  {
    Eigen::VectorXd x_bad = (Eigen::VectorXd(1) << 0.0).finished();
    Eigen::VectorXd y_bad = (Eigen::VectorXd(1) << 10.0).finished();
    Eigen::VectorXd dydx_bad = (Eigen::VectorXd(1) << 1.0).finished();
    EXPECT_THROW(stan::math::interpolate_1d_cubic_hermite_setup(x_bad, y_bad, dydx_bad),
                 std::domain_error);
  }

  // Non-finite x
  {
    Eigen::VectorXd x_nan = (Eigen::VectorXd(4) << 0.0, std::nan(""), 2.5, 4.0).finished();
    EXPECT_THROW(stan::math::interpolate_1d_cubic_hermite_setup(x_nan, kTestY, kTestDydx),
                 std::domain_error);
  }

  // Non-ordered x
  {
    Eigen::VectorXd x_bad = (Eigen::VectorXd(4) << 0.0, 1.0, 1.0, 4.0).finished();
    EXPECT_THROW(stan::math::interpolate_1d_cubic_hermite_setup(x_bad, kTestY, kTestDydx),
                 std::domain_error);
  }

  // Non-finite y
  {
    Eigen::VectorXd y_nan = (Eigen::VectorXd(4) << 10.0, std::nan(""), 9.0, 11.0).finished();
    EXPECT_THROW(stan::math::interpolate_1d_cubic_hermite_setup(kTestX, y_nan, kTestDydx),
                 std::domain_error);
  }

  // Non-finite derivatives
  {
    Eigen::VectorXd dydx_inf = (Eigen::VectorXd(4) << 1.0, std::numeric_limits<double>::infinity(), 
                                                      -0.5, 1.5).finished();
    EXPECT_THROW(stan::math::interpolate_1d_cubic_hermite_setup(kTestX, kTestY, dydx_inf),
                 std::domain_error);
  }
}

// ============================================================================
// SECTION 5: Exception Handling Tests for Evaluation Function
// ============================================================================

/**
 * TEST: interpolate1DCubicHermiteEvalThrows
 * 
 * Verifies that interpolate_1d_cubic_spline throws std::domain_error
 * for invalid evaluation points.
 */
TEST(MathPrimFun, interpolate1DCubicHermiteEvalThrows) {
  auto W = stan::math::interpolate_1d_cubic_hermite_setup(kTestX, kTestY, kTestDydx);

  // Out of bounds (below)
  EXPECT_THROW(stan::math::interpolate_1d_cubic_spline(-1.0, W), std::domain_error);

  // Out of bounds (above)
  EXPECT_THROW(stan::math::interpolate_1d_cubic_spline(5.0, W), std::domain_error);

  // Non-finite: NaN
  EXPECT_THROW(stan::math::interpolate_1d_cubic_spline(std::nan(""), W), std::domain_error);

  // Non-finite: infinity
  EXPECT_THROW(stan::math::interpolate_1d_cubic_spline(std::numeric_limits<double>::infinity(), W),
               std::domain_error);
}

// ============================================================================
// SECTION 6: Edge Cases and Robustness Tests
// ============================================================================

/**
 * TEST: interpolate1DCubicHermiteEvalBoundaryConditions
 * 
 * Verifies that evaluation exactly at the boundaries works correctly.
 */
TEST(MathPrimFun, interpolate1DCubicHermiteEvalBoundaryConditions) {
  auto W = stan::math::interpolate_1d_cubic_hermite_setup(kTestX, kTestY, kTestDydx);
  auto boost_interp = create_boost_cubic_hermite_interpolator(kTestX, kTestY, kTestDydx);

  // At lower boundary
  double x_lower = kTestX(0);
  double expected_lower = boost_interp(x_lower);
  double computed_lower = stan::math::interpolate_1d_cubic_spline(x_lower, W);
  EXPECT_NEAR(expected_lower, computed_lower, 1e-10);

  // At upper boundary
  double x_upper = kTestX(kTestX.size() - 1);
  double expected_upper = boost_interp(x_upper);
  double computed_upper = stan::math::interpolate_1d_cubic_spline(x_upper, W);
  EXPECT_NEAR(expected_upper, computed_upper, 1e-10);
}

/**
 * TEST: interpolate1DCubicHermiteEvalVerySmallInterval
 * 
 * Tests robustness with 1e-6 small intervals between knots.
 * Extremely challenging: interval of 1e-6 requires very careful handling
 * of numerical precision in normalized parameter computation.
 */
TEST(MathPrimFun, interpolate1DCubicHermiteEvalVerySmallInterval) {
  Eigen::VectorXd x = (Eigen::VectorXd(4) << 0.0, 1e-6, 2e-6, 1.0).finished();
  Eigen::VectorXd y = (Eigen::VectorXd(4) << 0.0, 0.5, 1.0, 0.5).finished();
  Eigen::VectorXd dydx = (Eigen::VectorXd(4) << 1.0, 2.0, 1.0, 0.5).finished();

  auto W = stan::math::interpolate_1d_cubic_hermite_setup(x, y, dydx);
  auto boost_interp = create_boost_cubic_hermite_interpolator(x, y, dydx);

  // Test at all knot points
  for (Eigen::Index i = 0; i < x.size(); ++i) {
    double expected = boost_interp(x(i));
    double computed = stan::math::interpolate_1d_cubic_spline(x(i), W);
    EXPECT_NEAR(expected, computed, 1e-8) << "Failed at knot " << i;
  }

  // Test at intermediate point in first tiny interval [0, 1e-6]
  double x_mid1 = 5e-7;
  double expected_mid1 = boost_interp(x_mid1);
  double computed_mid1 = stan::math::interpolate_1d_cubic_spline(x_mid1, W);
  EXPECT_NEAR(expected_mid1, computed_mid1, 1e-8) << "Failed at x = 5e-7";

  // Test at intermediate point in second tiny interval [1e-6, 2e-6]
  double x_mid2 = 1.5e-6;
  double expected_mid2 = boost_interp(x_mid2);
  double computed_mid2 = stan::math::interpolate_1d_cubic_spline(x_mid2, W);
  EXPECT_NEAR(expected_mid2, computed_mid2, 1e-8) << "Failed at x = 1.5e-6";

  // Test near beginning of very tight intervals
  double x_near_start = 1e-8;  // Close to but not at x(0)
  double expected_near = boost_interp(x_near_start);
  double computed_near = stan::math::interpolate_1d_cubic_spline(x_near_start, W);
  EXPECT_NEAR(expected_near, computed_near, 1e-8) << "Failed at x = 1e-8";

  // Test at quarter point in first interval
  double x_quarter = 2.5e-7;
  double expected_quarter = boost_interp(x_quarter);
  double computed_quarter = stan::math::interpolate_1d_cubic_spline(x_quarter, W);
  EXPECT_NEAR(expected_quarter, computed_quarter, 1e-8) << "Failed at x = 2.5e-7";
}

/**
 * TEST: interpolate1DCubicHermiteEvalZeroDerivatives
 * 
 * Tests behavior when all derivatives are zero (should reduce to cubic interpolation
 * with zero slopes at all knots).
 */
TEST(MathPrimFun, interpolate1DCubicHermiteEvalZeroDerivatives) {
  Eigen::VectorXd x = (Eigen::VectorXd(3) << 0.0, 1.0, 2.0).finished();
  Eigen::VectorXd y = (Eigen::VectorXd(3) << 10.0, 20.0, 15.0).finished();
  Eigen::VectorXd dydx = (Eigen::VectorXd(3) << 0.0, 0.0, 0.0).finished();

  auto W = stan::math::interpolate_1d_cubic_hermite_setup(x, y, dydx);
  auto boost_interp = create_boost_cubic_hermite_interpolator(x, y, dydx);

  double x_val = 0.5;
  double expected = boost_interp(x_val);
  double computed = stan::math::interpolate_1d_cubic_spline(x_val, W);
  EXPECT_NEAR(expected, computed, 1e-10);
}

/**
 * TEST: interpolate1DCubicHermiteEvalMixedSignDerivatives
 * 
 * Tests with mixed positive and negative derivatives.
 */
TEST(MathPrimFun, interpolate1DCubicHermiteEvalMixedSignDerivatives) {
  Eigen::VectorXd x = (Eigen::VectorXd(4) << 0.0, 1.0, 2.0, 3.0).finished();
  Eigen::VectorXd y = (Eigen::VectorXd(4) << 5.0, 10.0, 8.0, 12.0).finished();
  Eigen::VectorXd dydx = (Eigen::VectorXd(4) << 5.0, -2.0, 3.0, -1.0).finished();

  auto W = stan::math::interpolate_1d_cubic_hermite_setup(x, y, dydx);
  auto boost_interp = create_boost_cubic_hermite_interpolator(x, y, dydx);

  std::vector<double> test_points = {0.0, 0.33, 0.67, 1.0, 1.5, 2.0, 2.5, 3.0};
  
  for (double x_val : test_points) {
    double expected = boost_interp(x_val);
    double computed = stan::math::interpolate_1d_cubic_spline(x_val, W);
    EXPECT_NEAR(expected, computed, 1e-10) << "Failed at x = " << x_val;
  }
}

// ============================================================================
// SECTION 7: PCHIP Interpolation Tests
// ============================================================================

/**
 * Helper function to create a Boost PCHIP interpolator from Eigen vectors.
 */
auto create_boost_pchip_interpolator(const Eigen::VectorXd& x, const Eigen::VectorXd& y) {
  std::vector<double> x_vec(x.data(), x.data() + x.size());
  std::vector<double> y_vec(y.data(), y.data() + y.size());
  return boost::math::interpolators::pchip<std::vector<double>>(std::move(x_vec), std::move(y_vec));
}

/**
 * TEST: interpolate1DPchipEval
 * 
 * Verifies that interpolate_1d_pchip_setup correctly sets up PCHIP interpolation
 * and produces results matching Boost's PCHIP interpolator.
 */
TEST(MathPrimFun, interpolate1DPchipEval) {
  auto W = stan::math::interpolate_1d_pchip_setup(kTestX, kTestY);
  auto boost_interp = create_boost_pchip_interpolator(kTestX, kTestY);

  // Test at knot points
  for (Eigen::Index i = 0; i < kTestX.size(); ++i) {
    double expected = boost_interp(kTestX(i));
    double computed = stan::math::interpolate_1d_cubic_spline(kTestX(i), W);
    EXPECT_NEAR(expected, computed, 1e-10) << "Failed at knot " << i;
  }

  // Test at intermediate points
  std::vector<double> test_points = {0.5, 1.5, 2.0, 3.0};
  for (double x : test_points) {
    double expected = boost_interp(x);
    double computed = stan::math::interpolate_1d_cubic_spline(x, W);
    EXPECT_NEAR(expected, computed, 1e-10) << "Failed at x = " << x;
  }
}

// ============================================================================
// SECTION 8: MAKIMA Interpolation Tests
// ============================================================================

/**
 * Helper function to create a Boost MAKIMA interpolator from Eigen vectors.
 */
auto create_boost_makima_interpolator(const Eigen::VectorXd& x, const Eigen::VectorXd& y) {
  std::vector<double> x_vec(x.data(), x.data() + x.size());
  std::vector<double> y_vec(y.data(), y.data() + y.size());
  return boost::math::interpolators::makima<std::vector<double>>(std::move(x_vec), std::move(y_vec));
}

/**
 * TEST: interpolate1DMakimaEval
 *
 * Verifies that interpolate_1d_makima_setup correctly sets up MAKIMA
 * interpolation and produces results matching Boost's MAKIMA interpolator
 * at knot and intermediate points.
 */
TEST(MathPrimFun, interpolate1DMakimaEval) {
  auto W = stan::math::interpolate_1d_makima_setup(kTestX, kTestY);
  auto boost_interp = create_boost_makima_interpolator(kTestX, kTestY);

  // Test at knot points
  for (Eigen::Index i = 0; i < kTestX.size(); ++i) {
    double expected = boost_interp(kTestX(i));
    double computed = stan::math::interpolate_1d_cubic_spline(kTestX(i), W);
    EXPECT_NEAR(expected, computed, 1e-10) << "Failed at knot " << i;
  }

  // Test at intermediate points
  std::vector<double> test_points = {0.5, 1.5, 2.0, 3.0};
  for (double x : test_points) {
    double expected = boost_interp(x);
    double computed = stan::math::interpolate_1d_cubic_spline(x, W);
    EXPECT_NEAR(expected, computed, 1e-10) << "Failed at x = " << x;
  }
}


