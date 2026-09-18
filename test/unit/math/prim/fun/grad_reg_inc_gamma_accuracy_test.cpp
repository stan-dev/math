#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <string>

/**
 * Fixed-reference accuracy tests for the incomplete gamma gradient roots.
 *
 * Every expected value comes from mpmath at 80 digits, through two
 * independent routes that agree to at least 13 digits: numerical
 * differentiation of the regularized function, and an explicit integral
 * identity evaluated in density form. A third route, the Gautschi series,
 * confirms the deep-tail values. The generator is
 * `inc_gamma_beta_work/ref_mpmath.py` in the audit workspace.
 *
 * The references reproduce two constants that were already in
 * grad_reg_inc_gamma_test.cpp, which is the check that the two agree:
 * (a, z) = (1.1, 0.2) gives 0.31416364892410886, and (9, 10) gives
 * 0.120855166827777.
 *
 * These tests fail on develop at 5252d51d47. They exist because
 * `expect_ad` compares against finite differences and cannot see an error
 * of this size, and because the existing suites assert absolute tolerances
 * on a grid that stops short of every failing region.
 */

namespace {

/** Relative comparison, needed because the values span 1e-1 to 1e-148. */
void expect_rel(double expected, double actual, double tol,
                const std::string& label) {
  ASSERT_TRUE(std::isfinite(actual))
      << label << ": result is not finite, got " << actual;
  const double err = std::fabs(actual - expected) / std::fabs(expected);
  EXPECT_LT(err, tol) << label << ": relative error " << err << ", expected "
                      << expected << ", got " << actual;
}

double upper(double a, double z) {
  return stan::math::grad_reg_inc_gamma(a, z, stan::math::tgamma(a),
                                        stan::math::digamma(a));
}

}  // namespace

/**
 * D1. `z >= a && z >= 8` selects a fixed 10-term asymptotic expansion that
 * needs z >> a. On develop the error is 4.4e-2 at a = 20, 5.8e-1 at
 * a = 100 and 7.2e-1 at a = 171.
 */
TEST(MathPrimScalFun, grad_reg_inc_gamma_diagonal_z_equals_a) {
  expect_rel(8.95792055391225334e-02, upper(20.0, 20.0), 1e-10, "a=z=20");
  expect_rel(3.99274978577860498e-02, upper(100.0, 100.0), 1e-10, "a=z=100");
  expect_rel(3.05227525380366363e-02, upper(171.0, 171.0), 1e-10, "a=z=171");
}

/**
 * D2. Above a = 171.62, tgamma(a) is inf, so the correction term is
 * inf/inf. On develop every one of these returns NaN.
 */
TEST(MathPrimScalFun, grad_reg_inc_gamma_above_tgamma_overflow) {
  expect_rel(2.97491740225716214e-02, upper(180.0, 180.0), 1e-10, "a=z=180");
  expect_rel(2.30393429308686434e-02, upper(300.0, 300.0), 1e-10, "a=z=300");
  expect_rel(1.78442151466593514e-02, upper(500.0, 500.0), 1e-10, "a=z=500");
}

/**
 * D3. `grad_reg_lower_inc_gamma` runs its Gautschi branch for every z > 42
 * with a >= 12, and there `emz * (log_z * sum_a - sum_b)` cancels
 * completely. On develop it returns about 1e-14 of noise, with arbitrary
 * sign, at all four points.
 */
TEST(MathPrimScalFun, grad_reg_lower_inc_gamma_far_upper_tail) {
  expect_rel(-8.28941360740501778e-22,
             stan::math::grad_reg_lower_inc_gamma(50.0, 150.0), 1e-10,
             "a=50 z=150");
  expect_rel(-2.79418207415687613e-54,
             stan::math::grad_reg_lower_inc_gamma(50.0, 250.0), 1e-10,
             "a=50 z=250");
  expect_rel(-5.33817508735222058e-148,
             stan::math::grad_reg_lower_inc_gamma(50.0, 500.0), 1e-10,
             "a=50 z=500");
  expect_rel(-4.74445289261568464e-69,
             stan::math::grad_reg_lower_inc_gamma(170.0, 510.0), 1e-10,
             "a=170 z=510");
}

/**
 * The two functions compute the same quantity with opposite sign. On
 * develop they disagree by 72 % at a = z = 171. This test needs no external
 * reference at all.
 */
TEST(MathPrimScalFun, grad_reg_inc_gamma_agrees_with_lower) {
  for (double a : {5.0, 20.0, 50.0, 100.0, 171.0}) {
    const double up = upper(a, a);
    const double lo = stan::math::grad_reg_lower_inc_gamma(a, a);
    expect_rel(-up, lo, 1e-10, "a=z=" + std::to_string(a));
  }
}

/** Points that already worked, to show the change is not a regression. */
TEST(MathPrimScalFun, grad_reg_inc_gamma_unchanged_points) {
  expect_rel(3.14163648924108863e-01, upper(1.1, 0.2), 1e-10, "a=1.1 z=0.2");
  expect_rel(1.20855166827777014e-01, upper(9.0, 10.0), 1e-10, "a=9 z=10");
  expect_rel(3.32046796768609819e-11, upper(2.5, 30.0), 1e-10, "a=2.5 z=30");
}
