#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>
#include <cmath>

TEST_F(AgradRev, mathMixScalFun_std_normal_lccdf) {
  auto f = [](const auto& y) { return stan::math::std_normal_lccdf(y); };

  stan::test::expect_ad(f, 50.0);
  stan::test::expect_ad(f, 20.0 * stan::math::SQRT_TWO);
  stan::test::expect_ad(f, 5.5);
  stan::test::expect_ad(f, 0.0);
  stan::test::expect_ad(f, -0.15);
  stan::test::expect_ad(f, -1.14);
  stan::test::expect_ad(f, -3.00);
  stan::test::expect_ad(f, -10.00);
  // P1: exp(x2) overflow, reached through the lccdf reflection
  stan::test::expect_ad(f, -37.6);
  stan::test::expect_ad(f, -40.0);
  // P2: erfc(|s|)^2 underflow, s in (-20, -19.2103)
  stan::test::expect_ad(f, 27.5);
  stan::test::expect_ad(f, 28.0);
  // P3: quadratic residual fit diverges from the linear Mills asymptote
  stan::test::expect_ad(f, 60.0);
  stan::test::expect_ad(f, 100.0);

  // third order autodiff tests can fail at borders of piecewise function
  // stan::test::ad_tolerances tols;
  // tols.grad_hessian_grad_hessian_ = 1e1;
  // stan::test::expect_ad(tols, f, 0.1 * stan::math::SQRT_TWO);
}

TEST_F(AgradRev, mathMixScalFun_std_normal_lccdf_branch_cutoffs) {
  auto f = [](const auto& y) { return stan::math::std_normal_lccdf(y); };

  /* Bracket every interior cutoff of the piecewise value and derivative,
   * given here in scaled_y units so they line up with the table in
   * prim/prob/std_normal_lcdf.hpp. lccdf(y) reflects to lcdf(-y), so y =
   * -scaled_y * sqrt(2). The offset stays well clear of
   * finite_diff_grad_hessian_auto's stencil, which spans about 1.2e-5 * max(1,
   * |y|), so neither side straddles the seam.
   */
  const double cutoffs[] = {2.9,  2.5,  2.1,  1.5,  0.8,   0.1,  0.0,
                            -2.1, -3.9, -4.0, -7.0, -17.0, -29.0};
  for (double cut : cutoffs) {
    stan::test::expect_ad(f, -(cut - 0.01) * stan::math::SQRT_TWO);
    stan::test::expect_ad(f, -(cut + 0.01) * stan::math::SQRT_TWO);
  }
}

TEST_F(AgradRev, mathMixScalFun_std_normal_lccdf_branch_accuracy) {
  /* Pins the accuracy of each piecewise Taylor branch of the derivative
   * against d/dy log Phi(y) = phi(y)/Phi(y) at 60 significant
   * digits, by the mpmath recipe committed in
   * Phi(y) = erfc(-y/sqrt(2))/2, phi(y) = exp(-y^2/2)/sqrt(2 pi), rounded to
   * double. The tolerance on each row is that
   * branch's measured worst in-range relative error plus a small margin, so a
   * regression in any one branch names itself instead of hiding under
   * expect_ad's blanket 1e-4 relative gradient tolerance -- the (0.8, 1.5]
   * branch has only 1.6x margin there. These rows are the enforced half of the
   * cutoff table in prim/prob/std_normal_lcdf.hpp.
   */
  struct branch_ref {
    double y;
    double d1;
    double tol;
  };
  const branch_ref cases[] = {
      // (2.5, 2.9] of scaled_y, Taylor centre 2.7, measured worst 3.27e-05
      {-3.5496760415564683, -0.0007326476540693806, 5e-5},
      {-3.818376618407357, -0.00027222779390121885, 5e-5},
      {-4.1012193308819755, -8.881828792625139e-05, 5e-5},
      // (2.1, 2.5] of scaled_y, Taylor centre 2.3, measured worst 2.80e-05
      {-2.9839906166072305, -0.004655923761682003, 4e-5},
      {-3.2526911934581184, -0.00201252166907532, 4e-5},
      {-3.5355339059327378, -0.0007702965121768906, 4e-5},
      // (1.5, 2.1] of scaled_y, Taylor centre 1.85, measured worst 2.30e-05
      {-2.135462479183374, -0.04148009614764338, 3e-5},
      {-2.5455844122715714, -0.015709826784999336, 3e-5},
      {-2.9698484809835, -0.004856449376114489, 3e-5},
      // (0.8, 1.5] of scaled_y, Taylor centre 1.15, measured worst 6.09e-05
      {-1.145512985522207, -0.236841175835376, 9e-5},
      {-1.6263455967290592, -0.1121292480767062, 9e-5},
      {-2.121320343559643, -0.042773100995777136, 9e-5},
      // (0.1, 0.8] of scaled_y, Taylor centre 0.45, measured worst 7.61e-06
      {-0.15556349186104046, -0.7015595130271106, 1e-5},
      {-0.6363961030678928, -0.4416330793820557, 1e-5},
      {-1.1313708498984762, -0.24150063210766093, 1e-5},
  };

  for (const branch_ref& c : cases) {
    stan::math::var y = c.y;
    stan::math::var lp = stan::math::std_normal_lccdf(y);
    lp.grad();
    EXPECT_LT(std::fabs(y.adj() / c.d1 - 1.0), c.tol)
        << "derivative branch drifted at y = " << c.y;
    stan::math::set_zero_all_adjoints();
  }
}

namespace {

/** d2/dy2 log(1 - Phi(y)) */
double d2_std_normal_lccdf(double y) {
  stan::math::fvar<stan::math::fvar<double>> yv;
  yv.val_.val_ = y;
  yv.val_.d_ = 1.0;
  yv.d_.val_ = 1.0;
  return stan::math::std_normal_lccdf(yv).d_.d_;
}

/** d3/dy3 log(1 - Phi(y)) */
double d3_std_normal_lccdf(double y) {
  stan::math::fvar<stan::math::fvar<stan::math::fvar<double>>> yv;
  yv.val_.val_.val_ = y;
  yv.val_.val_.d_ = 1.0;
  yv.val_.d_.val_ = 1.0;
  yv.d_.val_.val_ = 1.0;
  return stan::math::std_normal_lccdf(yv).d_.d_.d_;
}

}  // namespace

/**
 * Lower-tail second and third derivatives. lccdf reflects to lcdf, so these
 * reuse the references in mix/prob/normal_cdf_log_test.cpp at -y, negated for
 * odd order. See that file for why they are asserted against a high-precision
 * reference rather than handed to expect_ad.
 */
TEST_F(AgradRev, mathMixScalFun_std_normal_lccdf_tail_derivatives) {
  struct tail_ref {
    double y;
    double expected;
  };
  const tail_ref d2_cases[] = {{-26.0, -1.67628759256344147e-146},
                               {-27.0, -5.39430099975435425e-158},
                               {-30.0, -4.42093840463564223e-195},
                               {-33.0, -4.42668490225313907e-236},
                               {-37.0, -7.84402424064104056e-297}};
  for (const tail_ref& c : d2_cases) {
    const double got = d2_std_normal_lccdf(c.y);
    ASSERT_FALSE(std::isnan(got)) << "NaN d2 at y = " << c.y;
    EXPECT_NE(0.0, got) << "d2 collapsed to zero at y = " << c.y;
    EXPECT_LT(std::fabs(got / c.expected - 1.0), 1e-6) << "d2 at y = " << c.y;
  }

  const tail_ref d3_cases[] = {{-20.0, -2.20285839650174548e-85},
                               {-22.0, -1.53317799001460164e-103},
                               {-25.0, -4.77605215552570066e-134},
                               {-30.0, -1.32480787525581414e-193}};
  for (const tail_ref& c : d3_cases) {
    const double got = d3_std_normal_lccdf(c.y);
    ASSERT_FALSE(std::isnan(got)) << "NaN d3 at y = " << c.y;
    EXPECT_NE(0.0, got) << "d3 collapsed to zero at y = " << c.y;
    EXPECT_LT(std::fabs(got / c.expected - 1.0), 1e-6) << "d3 at y = " << c.y;
  }
}

/** Sweep for non-finite derivatives across the whole useful range. */
TEST_F(AgradRev, mathMixScalFun_std_normal_lccdf_derivatives_are_finite) {
  for (double y = -120.0; y <= 120.0; y += 0.5) {
    EXPECT_FALSE(std::isnan(d2_std_normal_lccdf(y))) << "d2 NaN at y = " << y;
    EXPECT_FALSE(std::isnan(d3_std_normal_lccdf(y))) << "d3 NaN at y = " << y;
  }
}

/**
 * Far upper tail gradient. lccdf(y) reflects to lcdf(-y), so the gradient is
 * the inverse Mills ratio at -y with the sign flipped.
 */
TEST_F(AgradRev, mathMixScalFun_std_normal_lccdf_far_tail_gradient) {
  struct far_ref {
    double y;
    double value;
    double grad;
  };
  const far_ref cases[]
      = {{60.0, -1.80501356068056725e+03, -6.00166574202411240e+01},
         {100.0, -5.00552420869420530e+03, -1.00009998000999261e+02},
         {300.0, -4.50066227321186598e+04, -3.00003333259263400e+02}};
  for (const far_ref& c : cases) {
    stan::math::var y = c.y;
    stan::math::var lp = stan::math::std_normal_lccdf(y);
    lp.grad();
    EXPECT_LT(std::fabs(lp.val() / c.value - 1.0), 1e-12)
        << "value at y = " << c.y;
    EXPECT_LT(std::fabs(y.adj() / c.grad - 1.0), 1e-10)
        << "gradient at y = " << c.y;
    stan::math::set_zero_all_adjoints();
  }

  for (double y = 42.0; y <= 300.0; y += 0.25) {
    stan::math::var yv = y;
    stan::math::var lp = stan::math::std_normal_lccdf(yv);
    lp.grad();
    const double s = -y / std::sqrt(2.0);
    const double inv = 1.0 / (s * s);
    const double mills
        = -2.0 * s
          / (1.0 + inv * (-0.5 + inv * (0.75 + inv * (-1.875 + inv * 6.5625))));
    // dnlcdf is d/ds; the partial is dnlcdf * INV_SQRT_TWO
    EXPECT_LT(std::fabs(yv.adj() / (-mills / std::sqrt(2.0)) - 1.0), 1e-9)
        << "gradient drifts from the asymptotic Mills ratio at y = " << y;
    stan::math::set_zero_all_adjoints();
  }
}
