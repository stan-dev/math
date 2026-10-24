#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>
#include <gtest/gtest.h>
#include <tuple>
#include <cmath>

namespace {

// Fixed knot data (double scalars) shared across Hermite / PCHIP / MAKIMA
const Eigen::VectorXd kXk
    = (Eigen::VectorXd(5) << -1.0, 0.0, 1.2, 2.3, 3.0).finished();
const Eigen::VectorXd kYk
    = (Eigen::VectorXd(5) << 2.0, 4.0, 3.5, 6.0, 5.0).finished();
const Eigen::VectorXd kDydxk
    = (Eigen::VectorXd(5) << 1.5, 0.2, -0.8, 1.1, -0.3).finished();

}  // namespace

// ============================================================================
// expect_ad tests — evaluation at mid-interval points only.
// At the knots the second derivative jumps, so Hessian checks would fail.
// ============================================================================

/**
 * TEST: interpolate1DCubicHermiteExpectAD
 *
 * Automatic differentiation of `interpolate_1d_cubic_spline` w.r.t. the
 * evaluation point x, using a Hermite setup (user-provided derivatives).
 * Test points lie strictly between knots.
 */
TEST(MathMixFun, interpolate1DCubicHermiteExpectAD) {
  auto W = stan::math::interpolate_1d_cubic_hermite_setup(kXk, kYk, kDydxk);
  auto f = [&W](const auto& x) {
    return stan::math::interpolate_1d_cubic_spline(x, W);
  };

  // Midpoints of each of the 4 intervals (strictly between knots)
  stan::test::expect_ad(f, (kXk(0) + kXk(1)) / 2.0);  // interval 0
  stan::test::expect_ad(f, (kXk(1) + kXk(2)) / 2.0);  // interval 1
  stan::test::expect_ad(f, (kXk(2) + kXk(3)) / 2.0);  // interval 2
  stan::test::expect_ad(f, (kXk(3) + kXk(4)) / 2.0);  // interval 3
}

/**
 * TEST: interpolate1DPCHIPExpectAD
 *
 * Automatic differentiation of `interpolate_1d_cubic_spline` w.r.t. the
 * evaluation point x, using a PCHIP setup.  Test points lie strictly
 * between knots.
 */
TEST(MathMixFun, interpolate1DPCHIPExpectAD) {
  auto W = stan::math::interpolate_1d_pchip_setup(kXk, kYk);
  auto f = [&W](const auto& x) {
    return stan::math::interpolate_1d_cubic_spline(x, W);
  };

  stan::test::expect_ad(f, (kXk(0) + kXk(1)) / 2.0);
  stan::test::expect_ad(f, (kXk(1) + kXk(2)) / 2.0);
  stan::test::expect_ad(f, (kXk(2) + kXk(3)) / 2.0);
  stan::test::expect_ad(f, (kXk(3) + kXk(4)) / 2.0);
}

/**
 * TEST: interpolate1DMAKIMAExpectAD
 *
 * Automatic differentiation of `interpolate_1d_cubic_spline` w.r.t. the
 * evaluation point x, using a MAKIMA setup.  Test points lie strictly
 * between knots.
 */
TEST(MathMixFun, interpolate1DMAKIMAExpectAD) {
  auto W = stan::math::interpolate_1d_makima_setup(kXk, kYk);
  auto f = [&W](const auto& x) {
    return stan::math::interpolate_1d_cubic_spline(x, W);
  };

  stan::test::expect_ad(f, (kXk(0) + kXk(1)) / 2.0);
  stan::test::expect_ad(f, (kXk(1) + kXk(2)) / 2.0);
  stan::test::expect_ad(f, (kXk(2) + kXk(3)) / 2.0);
  stan::test::expect_ad(f, (kXk(3) + kXk(4)) / 2.0);
}