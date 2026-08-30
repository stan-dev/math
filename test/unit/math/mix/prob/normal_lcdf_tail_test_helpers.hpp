#ifndef TEST_UNIT_MATH_MIX_PROB_NORMAL_LCDF_TAIL_TEST_HELPERS_HPP
#define TEST_UNIT_MATH_MIX_PROB_NORMAL_LCDF_TAIL_TEST_HELPERS_HPP

#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>
#include <gtest/gtest.h>
#include <cmath>

/**
 * Shared checks for the tails of normal_lcdf, std_normal_lcdf and their
 * reflected counterparts normal_lccdf and std_normal_lccdf.
 *
 * All four share one piecewise implementation, so they share one set of
 * reference values. Those are given here for `normal_lcdf(y, 0, 1)`. A
 * reflected function -- `lccdf(y) == lcdf(-y)` -- is exercised at `-y`, and its
 * odd-order derivatives flip sign while even-order ones do not. Pass
 * `orientation::lcdf` or `orientation::lccdf` to select.
 *
 * References were evaluated at 60 significant digits from
 * `Phi(y) = erfc(-y/sqrt(2))/2` and `phi(y) = exp(-y^2/2)/sqrt(2 pi)`, then
 * rounded to double. The cutoff tables they enforce are documented on
 * stan::math::normal_lcdf in prim/prob/normal_lcdf.hpp.
 */
namespace normal_lcdf_tail_test {

/** Sign relating the function under test to normal_lcdf. */
struct orientation {
  static constexpr double lcdf = 1.0;
  static constexpr double lccdf = -1.0;
};

/** d2/dy2 of `f` at `y`, pure forward. */
template <typename F>
double d2(const F& f, double y) {
  stan::math::fvar<stan::math::fvar<double>> yv;
  yv.val_.val_ = y;
  yv.val_.d_ = 1.0;
  yv.d_.val_ = 1.0;
  return f(yv).d_.d_;
}

/**
 * d3/dy3 of `f` at `y`, pure forward.
 *
 * `fvar<fvar<fvar<double>>>` is deliberate. test/unit/math/test_ad.hpp defines
 * `d_t, v_t, fd_t, ffd_t, fv_t, ffv_t` and stops at `fvar<fvar<var>>`, so
 * expect_ad never instantiates this type at any input.
 */
template <typename F>
double d3(const F& f, double y) {
  stan::math::fvar<stan::math::fvar<stan::math::fvar<double>>> yv;
  yv.val_.val_.val_ = y;
  yv.val_.val_.d_ = 1.0;
  yv.val_.d_.val_ = 1.0;
  yv.d_.val_.val_ = 1.0;
  return f(yv).d_.d_.d_;
}

/** Value and reverse-mode gradient of `f` at `y`. */
template <typename F>
void value_and_grad(const F& f, double y, double* value, double* grad) {
  stan::math::var yv = y;
  stan::math::var lp = f(yv);
  lp.grad();
  *value = lp.val();
  *grad = yv.adj();
  stan::math::set_zero_all_adjoints();
}

/**
 * Five-term DLMF 7.12.1 asymptotic for d/dy log Phi(y), valid in the far lower
 * tail. Independent of the four-term form the implementation uses.
 */
inline double mills_dlogphi_dy(double y) {
  const double s = y / std::sqrt(2.0);
  const double inv = 1.0 / (s * s);
  const double series
      = 1.0 + inv * (-0.5 + inv * (0.75 + inv * (-1.875 + inv * 6.5625)));
  return -2.0 * s / series / std::sqrt(2.0);
}

/**
 * expect_ad on both sides of every interior cutoff of the piecewise value and
 * derivative, given in scaled units so they line up with the table in
 * prim/prob/normal_lcdf.hpp.
 *
 * The 0.01 offset keeps clear of `finite_diff_grad_hessian_auto`'s stencil,
 * which spans about `1.2e-5 * max(1, |y|)`, so neither side straddles a seam.
 *
 * This pins that every branch is executed and is accurate where used. It
 * cannot pin the cutoff location: adjacent branches agree to well within 1e-4
 * for about 0.1 either side of a seam, by construction. Mutation-checked --
 * re-centring the 1.85 series on 1.80 fails, shifting the 2.5 boundary to 2.4
 * does not, 2.75 does.
 */
template <typename F>
void expect_ad_across_cutoffs(const F& f, double dir) {
  const double cutoffs[] = {2.9,  2.5,  2.1,  1.5,  0.8,   0.1,  0.0,
                            -2.1, -3.9, -4.0, -7.0, -17.0, -29.0};
  for (double cut : cutoffs) {
    stan::test::expect_ad(f, dir * (cut - 0.01) * stan::math::SQRT_TWO);
    stan::test::expect_ad(f, dir * (cut + 0.01) * stan::math::SQRT_TWO);
  }
}

/**
 * expect_ad at the inputs that exposed each of the three tail defects: the
 * `exp(x2)` overflow above `scaled_diff = 26.6`, the `erfc^2` underflow window
 * `scaled_diff in (-20, -19.2103)`, and the far-tail residual fit that
 * diverges from the linear Mills asymptote.
 */
template <typename F>
void expect_ad_at_defect_inputs(const F& f, double dir) {
  const double inputs[] = {37.6, 40.0, 50.0, -27.5, -28.0, -60.0, -100.0};
  for (double y : inputs) {
    stan::test::expect_ad(f, dir * y);
  }
}

/**
 * Per-branch derivative accuracy, each row at that branch's own measured worst
 * in-range error plus a small margin. Per-branch rather than blanket because
 * expect_ad's 1e-4 `gradient_grad_` leaves the (0.8, 1.5] branch only 1.6x.
 */
template <typename F>
void expect_branch_accuracy(const F& f, double dir) {
  struct branch_ref {
    double y;
    double d1;
    double tol;
  };
  const branch_ref cases[]
      = {// (2.5, 2.9], Taylor centre 2.7, measured worst 3.27e-05
         {3.5496760415564683, 0.0007326476540693806, 5e-5},
         {3.818376618407357, 0.00027222779390121885, 5e-5},
         {4.1012193308819755, 8.881828792625139e-05, 5e-5},
         // (2.1, 2.5], Taylor centre 2.3, measured worst 2.80e-05
         {2.9839906166072305, 0.004655923761682003, 4e-5},
         {3.2526911934581184, 0.00201252166907532, 4e-5},
         {3.5355339059327378, 0.0007702965121768906, 4e-5},
         // (1.5, 2.1], Taylor centre 1.85, measured worst 2.30e-05
         {2.135462479183374, 0.04148009614764338, 3e-5},
         {2.5455844122715714, 0.015709826784999336, 3e-5},
         {2.9698484809835, 0.004856449376114489, 3e-5},
         // (0.8, 1.5], Taylor centre 1.15, measured worst 6.09e-05
         {1.145512985522207, 0.236841175835376, 9e-5},
         {1.6263455967290592, 0.1121292480767062, 9e-5},
         {2.121320343559643, 0.042773100995777136, 9e-5},
         // (0.1, 0.8], Taylor centre 0.45, measured worst 7.61e-06
         {0.15556349186104046, 0.7015595130271106, 1e-5},
         {0.6363961030678928, 0.4416330793820557, 1e-5},
         {1.1313708498984762, 0.24150063210766093, 1e-5}};

  for (const branch_ref& c : cases) {
    const double y = dir * c.y;
    double value = 0;
    double grad = 0;
    value_and_grad(f, y, &value, &grad);
    EXPECT_LT(std::fabs(grad / (dir * c.d1) - 1.0), c.tol)
        << "derivative branch drifted at y = " << y;
  }
}

/**
 * Second and third derivatives in the tail against 60-digit references.
 *
 * Asserted directly rather than through expect_ad because expect_ad cannot see
 * either failure. The order-2 rows were silently zero before the `exp(-x2)`
 * rearrangement, and a true value of ~1e-158 is invisible to a
 * finite-difference comparison whose relative tolerance floors its denominator
 * at 1. The order-3 rows use a type expect_ad never instantiates, and fail
 * earlier than anything it does build -- silently 0 at 20, NaN from 22,
 * against 37.6 for order 2.
 */
template <typename F>
void expect_tail_derivatives(const F& f, double dir) {
  struct tail_ref {
    double y;
    double expected;
  };
  const tail_ref d2_cases[] = {{26.0, -1.67628759256344147e-146},
                               {27.0, -5.39430099975435425e-158},
                               {30.0, -4.42093840463564223e-195},
                               {33.0, -4.42668490225313907e-236},
                               {37.0, -7.84402424064104056e-297}};
  for (const tail_ref& c : d2_cases) {
    const double y = dir * c.y;
    const double got = d2(f, y);
    ASSERT_FALSE(std::isnan(got)) << "NaN d2 at y = " << y;
    EXPECT_NE(0.0, got) << "d2 collapsed to zero at y = " << y;
    EXPECT_LT(std::fabs(got / c.expected - 1.0), 1e-6) << "d2 at y = " << y;
  }

  const tail_ref d3_cases[] = {{20.0, 2.20285839650174548e-85},
                               {22.0, 1.53317799001460164e-103},
                               {25.0, 4.77605215552570066e-134},
                               {30.0, 1.32480787525581414e-193}};
  for (const tail_ref& c : d3_cases) {
    const double y = dir * c.y;
    const double got = d3(f, y);
    ASSERT_FALSE(std::isnan(got)) << "NaN d3 at y = " << y;
    EXPECT_NE(0.0, got) << "d3 collapsed to zero at y = " << y;
    EXPECT_LT(std::fabs(got / (dir * c.expected) - 1.0), 1e-6)
        << "d3 at y = " << y;
  }
}

/**
 * Sweep for non-finite derivatives across the whole useful range. Cheap at 0.5
 * spacing for two direct fvar evaluations, where expect_ad -- six AD modes plus
 * finite differencing per point -- would not be. This single loop would have
 * caught every tail defect in this function's history at once.
 */
template <typename F>
void expect_derivatives_finite(const F& f) {
  for (double y = -120.0; y <= 120.0; y += 0.5) {
    EXPECT_FALSE(std::isnan(d2(f, y))) << "d2 NaN at y = " << y;
    EXPECT_FALSE(std::isnan(d3(f, y))) << "d3 NaN at y = " << y;
  }
}

/**
 * Far-tail value and gradient against 60-digit references, then a sweep
 * against the five-term Mills asymptotic over the region the implementation
 * handles with its own four-term form.
 */
template <typename F>
void expect_far_tail_gradient(const F& f, double dir) {
  struct far_ref {
    double y;
    double value;
    double grad;
  };
  const far_ref cases[]
      = {{-60.0, -1.80501356068056725e+03, 6.00166574202411240e+01},
         {-100.0, -5.00552420869420530e+03, 1.00009998000999261e+02},
         {-300.0, -4.50066227321186598e+04, 3.00003333259263400e+02}};
  for (const far_ref& c : cases) {
    const double y = dir * c.y;
    double value = 0;
    double grad = 0;
    value_and_grad(f, y, &value, &grad);
    EXPECT_LT(std::fabs(value / c.value - 1.0), 1e-12) << "value at y = " << y;
    EXPECT_LT(std::fabs(grad / (dir * c.grad) - 1.0), 1e-10)
        << "gradient at y = " << y;
  }

  for (double ly = -300.0; ly <= -42.0; ly += 0.25) {
    const double y = dir * ly;
    double value = 0;
    double grad = 0;
    value_and_grad(f, y, &value, &grad);
    EXPECT_LT(std::fabs(grad / (dir * mills_dlogphi_dy(ly)) - 1.0), 1e-9)
        << "gradient drifts from the asymptotic Mills ratio at y = " << y;
  }
}

}  // namespace normal_lcdf_tail_test

#endif
