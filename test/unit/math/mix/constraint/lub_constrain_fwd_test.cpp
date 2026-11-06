#include <stan/math/mix.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <vector>

namespace lub_constrain_fwd_test {

// For lb = -2, ub = 3: x, lub_constrain(x, lb, ub), its derivatives with
// respect to lb, ub and x, its second derivative with respect to x, and the
// first and second derivatives of the log Jacobian with respect to x.
// mpmath, 50 digits, rounded to double.
struct ref {
  double x;
  double val;
  double d_lb;
  double d_ub;
  double d_x;
  double d_xx;
  double lp_x;
  double lp_xx;
};

// clang-format off
// NOLINTBEGIN(whitespace/line_length)
inline std::vector<ref> refs() {
  return {
      {0, 0.5, 0.5, 0.5, 1.25, 0, 0, -0.5},
      {1, 1.6552928931500244, 0.2689414213699951, 0.7310585786300049, 0.9830596662074093, -0.45428873836474204, -0.46211715726000974, -0.3932238664829637},
      {-1, -0.6552928931500244, 0.7310585786300049, 0.2689414213699951, 0.9830596662074093, 0.45428873836474204, 0.46211715726000974, -0.3932238664829637},
      {20, 2.9999999896942318, 2.0611536181902037e-09, 0.9999999979388464, 1.0305768069709247e-08, -1.0305768027225704e-08, -0.9999999958776927, -4.122307227883698e-09},
      {-20, -1.999999989694232, 0.9999999979388464, 2.0611536181902037e-09, 1.0305768069709247e-08, 1.0305768027225704e-08, 0.9999999958776927, -4.122307227883698e-09},
      {37, 2.9999999999999996, 8.533047625744065e-17, 0.9999999999999999, 4.2665238128720323e-16, -4.2665238128720313e-16, -0.9999999999999998, -1.706609525148813e-16},
      {-37, -1.9999999999999996, 0.9999999999999999, 8.533047625744065e-17, 4.2665238128720323e-16, 4.2665238128720313e-16, 0.9999999999999998, -1.706609525148813e-16},
      {40, 3.0, 4.248354255291589e-18, 1.0, 2.1241771276457944e-17, -2.1241771276457944e-17, -1.0, -8.496708510583178e-18},
      {-40, -2.0, 1.0, 4.248354255291589e-18, 2.1241771276457944e-17, 2.1241771276457944e-17, 1.0, -8.496708510583178e-18},
      {100, 3.0, 3.720075976020836e-44, 1.0, 1.860037988010418e-43, -1.860037988010418e-43, -1.0, -7.440151952041672e-44},
      {-100, -2.0, 1.0, 3.720075976020836e-44, 1.860037988010418e-43, 1.860037988010418e-43, 1.0, -7.440151952041672e-44},
      {700, 3.0, 9.85967654375977e-305, 1.0, 4.9298382718798856e-304, -4.9298382718798856e-304, -1.0, -1.971935308751954e-304},
      {-700, -2.0, 1.0, 9.85967654375977e-305, 4.9298382718798856e-304, 4.9298382718798856e-304, 1.0, -1.971935308751954e-304}};
}
// NOLINTEND
// clang-format on

// Relative check, absolute below the smallest normal double.
inline void expect_rel(double expected, double actual, double x,
                       const char* what) {
  EXPECT_NEAR(
      expected, actual,
      1e-13
          * std::fmax(std::fabs(expected), std::numeric_limits<double>::min()))
      << what << " at x = " << x;
}

}  // namespace lub_constrain_fwd_test

TEST(mathMixConstraint, lub_constrain_fwd_derivatives_tails) {
  using lub_constrain_fwd_test::expect_rel;
  using stan::math::fvar;
  using stan::math::lub_constrain;
  using stan::math::var;
  const double lb = -2.0;
  const double ub = 3.0;
  for (const auto& r : lub_constrain_fwd_test::refs()) {
    expect_rel(r.val, lub_constrain(r.x, lb, ub), r.x, "value");
    expect_rel(r.d_x, lub_constrain(fvar<double>(r.x, 1), lb, ub).d_, r.x,
               "d/dx");
    expect_rel(r.d_lb, lub_constrain(r.x, fvar<double>(lb, 1), ub).d_, r.x,
               "d/dlb");
    expect_rel(r.d_ub, lub_constrain(r.x, lb, fvar<double>(ub, 1)).d_, r.x,
               "d/dub");

    fvar<fvar<double>> x_ff(fvar<double>(r.x, 1), fvar<double>(1, 0));
    fvar<fvar<double>> y_ff = lub_constrain(x_ff, lb, ub);
    expect_rel(r.d_x, y_ff.d_.val_, r.x, "d/dx, fvar<fvar<double>>");
    expect_rel(r.d_xx, y_ff.d_.d_, r.x, "d2/dx2, fvar<fvar<double>>");

    fvar<var> x_fv(r.x, 1);
    fvar<var> y_fv = lub_constrain(x_fv, lb, ub);
    y_fv.d_.grad();
    expect_rel(r.d_xx, x_fv.val_.adj(), r.x, "d2/dx2, fvar<var>");
    stan::math::recover_memory();
  }
}

TEST(mathMixConstraint, lub_constrain_fwd_derivatives_tails_lp) {
  using lub_constrain_fwd_test::expect_rel;
  using stan::math::fvar;
  using stan::math::lub_constrain;
  using stan::math::var;
  const double lb = -2.0;
  const double ub = 3.0;
  for (const auto& r : lub_constrain_fwd_test::refs()) {
    double lp_d = 0;
    expect_rel(r.val, lub_constrain(r.x, lb, ub, lp_d), r.x, "value, lp");

    fvar<double> lp_lb = 0;
    expect_rel(r.d_lb, lub_constrain(r.x, fvar<double>(lb, 1), ub, lp_lb).d_,
               r.x, "d/dlb, lp");

    fvar<fvar<double>> x_ff(fvar<double>(r.x, 1), fvar<double>(1, 0));
    fvar<fvar<double>> lp_ff = 0;
    fvar<fvar<double>> y_ff = lub_constrain(x_ff, lb, ub, lp_ff);
    expect_rel(r.d_x, y_ff.d_.val_, r.x, "d/dx, lp");
    expect_rel(r.d_xx, y_ff.d_.d_, r.x, "d2/dx2, lp");
    expect_rel(r.lp_x, lp_ff.d_.val_, r.x, "log Jacobian d/dx");
    expect_rel(r.lp_xx, lp_ff.d_.d_, r.x,
               "log Jacobian d2/dx2, fvar<fvar<double>>");

    fvar<var> x_fv(r.x, 1);
    fvar<var> lp_fv = 0;
    lub_constrain(x_fv, lb, ub, lp_fv);
    lp_fv.d_.grad();
    expect_rel(r.lp_xx, x_fv.val_.adj(), r.x, "log Jacobian d2/dx2, fvar<var>");
    stan::math::recover_memory();
  }
}
