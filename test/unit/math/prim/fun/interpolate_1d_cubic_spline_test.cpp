#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <tuple>
#include <cmath>

namespace {
const Eigen::VectorXd kTestX = (Eigen::VectorXd(4) << 0.0, 1.0, 2.5, 4.0).finished();
const Eigen::VectorXd kTestY = (Eigen::VectorXd(4) << 10.0, 12.0, 9.0, 11.0).finished();
const Eigen::VectorXd kTestDydx = (Eigen::VectorXd(4) << 1.0, 0.5, -0.5, 1.5).finished();
}  // namespace

// ============================================================================
// SECTION 1: Setup Function Tests
// ============================================================================

TEST(MathPrimFun, interpolate1DCubicHermiteSetupStorage) {
  auto W = stan::math::interpolate_1d_cubic_hermite_setup(kTestX, kTestY, kTestDydx);
  const auto& xk = std::get<0>(W);
  const auto& yk = std::get<1>(W);
  const auto& dydxk = std::get<2>(W);
  EXPECT_EQ(kTestX.size(), xk.size());
  EXPECT_EQ(kTestY.size(), yk.size());
  EXPECT_EQ(kTestDydx.size(), dydxk.size());
  for (Eigen::Index i = 0; i < kTestX.size(); ++i)
    EXPECT_FLOAT_EQ(kTestX(i), xk(i));
  for (Eigen::Index i = 0; i < kTestY.size(); ++i)
    EXPECT_FLOAT_EQ(kTestY(i), yk(i));
  for (Eigen::Index i = 0; i < kTestDydx.size(); ++i)
    EXPECT_FLOAT_EQ(kTestDydx(i), dydxk(i));
}

// ============================================================================
// SECTION 2: Knot Value Tests
// ============================================================================

TEST(MathPrimFun, interpolate1DCubicHermiteEvalAtKnots) {
  auto W = stan::math::interpolate_1d_cubic_hermite_setup(kTestX, kTestY, kTestDydx);
  for (Eigen::Index i = 0; i < kTestX.size(); ++i) {
    double expected_y = kTestY(i);
    double computed_y = stan::math::interpolate_1d_cubic_spline(kTestX(i), W);
    EXPECT_NEAR(expected_y, computed_y, 1e-10) << "Failed at knot " << i;
  }
}

// ============================================================================
// SECTION 3: Exception Handling — Setup Function
// ============================================================================

TEST(MathPrimFun, interpolate1DCubicHermiteSetupThrowsInvalidArgument) {
  {
    Eigen::VectorXd x_bad = (Eigen::VectorXd(3) << 0.0, 1.0, 2.5).finished();
    Eigen::VectorXd y_bad = (Eigen::VectorXd(4) << 10.0, 12.0, 9.0, 11.0).finished();
    Eigen::VectorXd dydx_good = kTestDydx;
    EXPECT_THROW(stan::math::interpolate_1d_cubic_hermite_setup(x_bad, y_bad, dydx_good),
                 std::invalid_argument);
  }
  {
    Eigen::VectorXd x_good = kTestX;
    Eigen::VectorXd y_good = kTestY;
    Eigen::VectorXd dydx_bad = (Eigen::VectorXd(3) << 1.0, 0.5, -0.5).finished();
    EXPECT_THROW(stan::math::interpolate_1d_cubic_hermite_setup(x_good, y_good, dydx_bad),
                 std::invalid_argument);
  }
}

TEST(MathPrimFun, interpolate1DCubicHermiteSetupThrowsDomainErrors) {
  {
    Eigen::VectorXd x_bad = (Eigen::VectorXd(1) << 0.0).finished();
    Eigen::VectorXd y_bad = (Eigen::VectorXd(1) << 10.0).finished();
    Eigen::VectorXd dydx_bad = (Eigen::VectorXd(1) << 1.0).finished();
    EXPECT_THROW(stan::math::interpolate_1d_cubic_hermite_setup(x_bad, y_bad, dydx_bad),
                 std::domain_error);
  }
  {
    Eigen::VectorXd x_nan = (Eigen::VectorXd(4) << 0.0, std::nan(""), 2.5, 4.0).finished();
    EXPECT_THROW(stan::math::interpolate_1d_cubic_hermite_setup(x_nan, kTestY, kTestDydx),
                 std::domain_error);
  }
  {
    Eigen::VectorXd x_bad = (Eigen::VectorXd(4) << 0.0, 1.0, 1.0, 4.0).finished();
    EXPECT_THROW(stan::math::interpolate_1d_cubic_hermite_setup(x_bad, kTestY, kTestDydx),
                 std::domain_error);
  }
  {
    Eigen::VectorXd y_nan = (Eigen::VectorXd(4) << 10.0, std::nan(""), 9.0, 11.0).finished();
    EXPECT_THROW(stan::math::interpolate_1d_cubic_hermite_setup(kTestX, y_nan, kTestDydx),
                 std::domain_error);
  }
  {
    Eigen::VectorXd dydx_inf = (Eigen::VectorXd(4) << 1.0,
        std::numeric_limits<double>::infinity(), -0.5, 1.5).finished();
    EXPECT_THROW(stan::math::interpolate_1d_cubic_hermite_setup(kTestX, kTestY, dydx_inf),
                 std::domain_error);
  }
}

// ============================================================================
// SECTION 4: Exception Handling — Evaluation Function
// ============================================================================

TEST(MathPrimFun, interpolate1DCubicHermiteEvalThrows) {
  auto W = stan::math::interpolate_1d_cubic_hermite_setup(kTestX, kTestY, kTestDydx);
  EXPECT_THROW(stan::math::interpolate_1d_cubic_spline(-1.0, W), std::domain_error);
  EXPECT_THROW(stan::math::interpolate_1d_cubic_spline(5.0, W), std::domain_error);
  EXPECT_THROW(stan::math::interpolate_1d_cubic_spline(std::nan(""), W), std::domain_error);
  EXPECT_THROW(stan::math::interpolate_1d_cubic_spline(
      std::numeric_limits<double>::infinity(), W), std::domain_error);
}

// ============================================================================
// SECTION 5: Property Tests — Reproduction of a Cubic
// ============================================================================

/**
 * TEST: interpolate1DCubicReproducesCubic
 *
 * If the data and derivatives come from an exact cubic polynomial, the
 * cubic Hermite interpolant must reproduce that cubic exactly at every
 * point, in every interval.  This test needs no reference values — it
 * compares against the generating polynomial f(x) = a x³ + b x² + c x + d.
 */
TEST(MathPrimFun, interpolate1DCubicReproducesCubic) {
  const double a = 2.0, b = -3.0, c = 1.0, d = 5.0;
  auto f = [&](double x) { return a * x * x * x + b * x * x + c * x + d; };
  auto fp = [&](double x) { return 3 * a * x * x + 2 * b * x + c; };

  // Knots: 5 non-uniformly spaced points spanning [-1.0, 3.0] → 4 intervals
  const int N = 5;
  Eigen::VectorXd xk(N), yk(N), dydxk(N);
  xk << -1.0, 0.0, 1.2, 2.3, 3.0;
  for (int i = 0; i < N; ++i) {
    yk(i) = f(xk(i));
    dydxk(i) = fp(xk(i));
  }

  auto W = stan::math::interpolate_1d_cubic_hermite_setup(xk, yk, dydxk);

  // Densely sample every interval
  for (int iv = 0; iv < N - 1; ++iv) {
    for (int k = 0; k <= 10; ++k) {
      double t = k / 10.0;
      double x = xk(iv) + t * (xk(iv + 1) - xk(iv));
      double spline_val = stan::math::interpolate_1d_cubic_spline(x, W);
      double true_val = f(x);
      EXPECT_NEAR(true_val, spline_val, 3e-13)
          << "Cubic reproduction failed at x=" << x << " in interval " << iv;
    }
  }
}

// ============================================================================
// SECTION 6: Property Tests — Continuity of Value and Slope at Knots
// ============================================================================

namespace {
template <typename W_t>
double numerical_deriv(double x, const W_t& W, double h = 1e-6) {
  return (stan::math::interpolate_1d_cubic_spline(x + h, W)
          - stan::math::interpolate_1d_cubic_spline(x - h, W)) / (2 * h);
}
}  // namespace

/**
 * TEST: interpolate1DCubicContinuityValueAndSlope
 *
 * At each interior knot the interpolated value must equal the specified
 * knot value, and the numerical derivative approached from left and
 * right must both match the specified derivative.
 */
TEST(MathPrimFun, interpolate1DCubicContinuityValueAndSlope) {
  auto W = stan::math::interpolate_1d_cubic_hermite_setup(kTestX, kTestY, kTestDydx);

  const double eps = 1e-6;
  for (Eigen::Index i = 1; i < kTestX.size() - 1; ++i) {
    double xk_i = kTestX(i);

    // Value continuity: approaching from left and right, both one-sided
    // limits must converge to the knot value.  At distance `eps` the value
    // deviates by ~ |dydx(i)| * eps, so the tolerance scales with eps.
    double val_tol = 2 * std::abs(kTestDydx(i)) * eps + 1e-9;
    double val_left = stan::math::interpolate_1d_cubic_spline(xk_i - eps, W);
    double val_right = stan::math::interpolate_1d_cubic_spline(xk_i + eps, W);
    EXPECT_NEAR(kTestY(i), val_left, val_tol)
        << "Value continuity (left) at knot " << i;
    EXPECT_NEAR(kTestY(i), val_right, val_tol)
        << "Value continuity (right) at knot " << i;

    // Slope continuity: numerical derivative from left and right
    double slope_left = numerical_deriv(xk_i - eps, W);
    double slope_right = numerical_deriv(xk_i + eps, W);
    EXPECT_NEAR(kTestDydx(i), slope_left, 5e-5)
        << "Slope continuity (left) at knot " << i;
    EXPECT_NEAR(kTestDydx(i), slope_right, 5e-5)
        << "Slope continuity (right) at knot " << i;
  }
}

// ============================================================================
// SECTION 7: Property Tests — PCHIP Monotonicity
// ============================================================================

/**
 * TEST: interpolate1DPCHIPMonotone
 *
 * When the data are strictly monotone, the PCHIP interpolant must be
 * monotone (no spurious oscillations).  This test uses strictly
 * increasing data and verifies the interpolant is strictly increasing
 * in every interval.
 */
TEST(MathPrimFun, interpolate1DPCHIPMonotone) {
  const int N = 6;
  Eigen::VectorXd xk = (Eigen::VectorXd(N) << 0.0, 1.0, 2.5, 4.0, 6.0, 8.0).finished();
  // Strictly increasing: y = x²
  Eigen::VectorXd yk(N);
  for (int i = 0; i < N; ++i)
    yk(i) = xk(i) * xk(i);

  auto W = stan::math::interpolate_1d_pchip_setup(xk, yk);

  // Sample densely in every interval
  for (int iv = 0; iv < N - 1; ++iv) {
    double prev = -std::numeric_limits<double>::infinity();
    for (int k = 0; k <= 20; ++k) {
      double t = k / 20.0;
      double x = xk(iv) + t * (xk(iv + 1) - xk(iv));
      double val = stan::math::interpolate_1d_cubic_spline(x, W);
      EXPECT_GT(val, prev) << "PCHIP not monotone at x=" << x;
      prev = val;
    }
  }

  // Also test strictly decreasing data
  Eigen::VectorXd yk_dec(N);
  for (int i = 0; i < N; ++i)
    yk_dec(i) = -xk(i) * xk(i);

  auto W_dec = stan::math::interpolate_1d_pchip_setup(xk, yk_dec);

  for (int iv = 0; iv < N - 1; ++iv) {
    double prev = std::numeric_limits<double>::infinity();
    for (int k = 0; k <= 20; ++k) {
      double t = k / 20.0;
      double x = xk(iv) + t * (xk(iv + 1) - xk(iv));
      double val = stan::math::interpolate_1d_cubic_spline(x, W_dec);
      EXPECT_LT(val, prev) << "PCHIP not monotone (decreasing) at x=" << x;
      prev = val;
    }
  }
}

// ============================================================================
// SECTION 8: Property Tests — Invariance under Scaling
// ============================================================================

/**
 * TEST: interpolate1DCubicScaleInvariance
 *
 * Under an affine transformation of x and y the interpolant must
 * transform predictably:
 *   x' = α x + β,  y' = γ y + δ  →  spline'(α x + β) == γ·spline(x) + δ
 *   dydx' = (γ / α) · dydx
 */
TEST(MathPrimFun, interpolate1DCubicScaleInvariance) {
  const double alpha = 3.0, beta = -2.0, gamma = 0.5, delta = 7.0;

  // Scale the knot data
  Eigen::VectorXd xk_scaled = alpha * kTestX.array() + beta;
  Eigen::VectorXd yk_scaled = gamma * kTestY.array() + delta;
  Eigen::VectorXd dydx_scaled = (gamma / alpha) * kTestDydx.array();

  auto W_base = stan::math::interpolate_1d_cubic_hermite_setup(kTestX, kTestY, kTestDydx);
  auto W_scaled = stan::math::interpolate_1d_cubic_hermite_setup(xk_scaled, yk_scaled, dydx_scaled);

  // Check at dense points in every interval of the base data
  const int N = kTestX.size();
  for (int iv = 0; iv < N - 1; ++iv) {
    for (int k = 0; k <= 10; ++k) {
      double t = k / 10.0;
      double x_base = kTestX(iv) + t * (kTestX(iv + 1) - kTestX(iv));
      double x_scaled = alpha * x_base + beta;

      double val_base = stan::math::interpolate_1d_cubic_spline(x_base, W_base);
      double val_scaled = stan::math::interpolate_1d_cubic_spline(x_scaled, W_scaled);
      double expected = gamma * val_base + delta;

      EXPECT_NEAR(expected, val_scaled, 5e-13)
          << "Scale invariance failed at x_base=" << x_base;
    }
  }
}

// ============================================================================
// SECTION 9: Invalid Input Tests for PCHIP Setup
// ============================================================================

/**
 * TEST: interpolate1DPCHIPSetupThrows
 *
 * Verifies that interpolate_1d_pchip_setup throws appropriate exceptions
 * for invalid inputs (size mismatch, too few knots, non-finite, non-ordered).
 */
TEST(MathPrimFun, interpolate1DPCHIPSetupThrows) {
  // Size mismatch
  {
    Eigen::VectorXd x(4), y(3);
    x << 0.0, 1.0, 2.0, 3.0;
    y << 0.0, 1.0, 2.0;
    EXPECT_THROW(stan::math::interpolate_1d_pchip_setup(x, y),
                 std::invalid_argument);
  }

  // Fewer than 4 knots
  {
    Eigen::VectorXd x(3), y(3);
    x << 0.0, 1.0, 2.0;
    y << 0.0, 1.0, 4.0;
    EXPECT_THROW(stan::math::interpolate_1d_pchip_setup(x, y),
                 std::domain_error);
  }

  // Non-finite x
  {
    Eigen::VectorXd x(4), y(4);
    x << 0.0, std::nan(""), 2.0, 3.0;
    y << 0.0, 1.0, 4.0, 9.0;
    EXPECT_THROW(stan::math::interpolate_1d_pchip_setup(x, y),
                 std::domain_error);
  }

  // Non-finite y
  {
    Eigen::VectorXd x(4), y(4);
    x << 0.0, 1.0, 2.0, 3.0;
    y << 0.0, std::numeric_limits<double>::infinity(), 4.0, 9.0;
    EXPECT_THROW(stan::math::interpolate_1d_pchip_setup(x, y),
                 std::domain_error);
  }

  // Non-ordered x
  {
    Eigen::VectorXd x(4), y(4);
    x << 0.0, 3.0, 2.0, 4.0;
    y << 0.0, 9.0, 4.0, 16.0;
    EXPECT_THROW(stan::math::interpolate_1d_pchip_setup(x, y),
                 std::domain_error);
  }
}

// ============================================================================
// SECTION 10: Invalid Input Tests for MAKIMA Setup
// ============================================================================

/**
 * TEST: interpolate1DMAKIMASetupThrows
 *
 * Verifies that interpolate_1d_makima_setup throws appropriate exceptions
 * for invalid inputs (size mismatch, too few knots, non-finite, non-ordered).
 */
TEST(MathPrimFun, interpolate1DMAKIMASetupThrows) {
  // Size mismatch
  {
    Eigen::VectorXd x(4), y(3);
    x << 0.0, 1.0, 2.0, 3.0;
    y << 0.0, 1.0, 2.0;
    EXPECT_THROW(stan::math::interpolate_1d_makima_setup(x, y),
                 std::invalid_argument);
  }

  // Fewer than 4 knots
  {
    Eigen::VectorXd x(3), y(3);
    x << 0.0, 1.0, 2.0;
    y << 0.0, 1.0, 4.0;
    EXPECT_THROW(stan::math::interpolate_1d_makima_setup(x, y),
                 std::domain_error);
  }

  // Non-finite x
  {
    Eigen::VectorXd x(4), y(4);
    x << 0.0, std::nan(""), 2.0, 3.0;
    y << 0.0, 1.0, 4.0, 9.0;
    EXPECT_THROW(stan::math::interpolate_1d_makima_setup(x, y),
                 std::domain_error);
  }

  // Non-finite y
  {
    Eigen::VectorXd x(4), y(4);
    x << 0.0, 1.0, 2.0, 3.0;
    y << 0.0, std::numeric_limits<double>::infinity(), 4.0, 9.0;
    EXPECT_THROW(stan::math::interpolate_1d_makima_setup(x, y),
                 std::domain_error);
  }

  // Non-ordered x
  {
    Eigen::VectorXd x(4), y(4);
    x << 0.0, 3.0, 2.0, 4.0;
    y << 0.0, 9.0, 4.0, 16.0;
    EXPECT_THROW(stan::math::interpolate_1d_makima_setup(x, y),
                 std::domain_error);
  }
}
