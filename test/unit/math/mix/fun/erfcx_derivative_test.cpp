#include <stan/math/mix.hpp>
#include <gtest/gtest.h>
#include <cmath>

/**
 * Accuracy of the erfcx derivative in the tail.
 *
 * `d/dx erfcx(x) = 2 * x * erfcx(x) - 2 / sqrt(pi)`. For large `x` the two
 * terms both approach `2 / sqrt(pi)` while their difference decays like
 * `1 / (sqrt(pi) * x^2)`, so the subtraction loses the result. Measured
 * against a 50-digit reference, the difference form gives 6.1 ulp at x = 4
 * and 2.55e+11 ulp at x = 1e6.
 *
 * For `x >= 4` the same tail rational that gives the value also gives the
 * derivative with no subtraction at all, because the `2 / sqrt(pi)` cancels
 * analytically:
 *
 *   erfcx(x)      = (INV_SQRT_PI + u * C(u)) / x,   u = 1 / x^2
 *   d/dx erfcx(x) = 2 * u * C(u)
 *
 * The existing fvar and var suites compare against finite differences, whose
 * tolerance is far looser than this error, so they pass either way. These
 * rows are fixed references instead.
 *
 * References were evaluated at 60 significant digits from
 * `2 * x * exp(x^2) * erfc(x) - 2 / sqrt(pi)` and rounded to double.
 */

namespace {

struct deriv_ref {
  double x;
  double d1;
};

// The suggested form measures 0.1 to 31 ulp across these points. 1e-13
// relative leaves room for that and still rejects the difference form, which
// is wrong by 6.8e+03 ulp at x = 100 and more beyond.
const deriv_ref TAIL[] = {{4.0, -0.032383506095021455},
                          {6.0, -0.015060353489052321},
                          {10.0, -0.0055593122190608565},
                          {100.0, -5.6410497625993184e-05},
                          {1000.0, -5.641887372654967e-07},
                          {10000.0, -5.641895750849127e-09}};

// Below 4 the difference form has no cancellation and stays in use.
const deriv_ref INTERIOR[] = {{-1.0, -11.14633932862008},
                              {-5.0, -1440097986747.7388}};

constexpr double TOL = 1e-13;

template <typename F>
void expect_rows(const deriv_ref* rows, int n, const F& derivative,
                 const char* what) {
  for (int i = 0; i < n; ++i) {
    const double got = derivative(rows[i].x);
    EXPECT_LT(std::fabs(got / rows[i].d1 - 1.0), TOL)
        << what << " derivative lost precision at x = " << rows[i].x;
  }
}

double fwd_derivative(double x) {
  stan::math::fvar<double> xv(x, 1.0);
  return stan::math::erfcx(xv).d_;
}

double rev_derivative(double x) {
  stan::math::var xv = x;
  stan::math::var y = stan::math::erfcx(xv);
  y.grad();
  const double d = xv.adj();
  stan::math::set_zero_all_adjoints();
  return d;
}

}  // namespace

TEST(MathFunctions, erfcx_derivative_tail_fwd) {
  expect_rows(TAIL, 6, fwd_derivative, "fvar tail");
}

TEST(MathFunctions, erfcx_derivative_tail_rev) {
  expect_rows(TAIL, 6, rev_derivative, "var tail");
}

TEST(MathFunctions, erfcx_derivative_interior_fwd) {
  expect_rows(INTERIOR, 2, fwd_derivative, "fvar interior");
}

TEST(MathFunctions, erfcx_derivative_interior_rev) {
  expect_rows(INTERIOR, 2, rev_derivative, "var interior");
}

/**
 * The two forms must agree where they meet, so the branch at 4 introduces no
 * step in the derivative.
 */
TEST(MathFunctions, erfcx_derivative_no_step_at_four) {
  const double below = fwd_derivative(std::nextafter(4.0, 0.0));
  const double at = fwd_derivative(4.0);
  EXPECT_LT(std::fabs(below / at - 1.0), 1e-12);
}
