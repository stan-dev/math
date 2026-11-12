#ifndef TEST_UNIT_MATH_MIX_FUN_INV_LOGIT_TAIL_REFS_HPP
#define TEST_UNIT_MATH_MIX_FUN_INV_LOGIT_TAIL_REFS_HPP

#include <stan/math/mix.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <vector>

namespace inv_logit_tail_refs {

// inv_logit(x), inv_logit(-x), inv_logit'(x) and inv_logit''(x), computed
// with mpmath at 60 digits and rounded to double (0 where they underflow).
struct ref {
  double x;
  double s;
  double s_neg;
  double d1;
  double d2;
};

inline std::vector<ref> refs() {
  return {{-800, 0, 1.0, 0, 0},
          {-40, 4.248354255291589e-18, 1.0, 4.248354255291589e-18,
           4.248354255291589e-18},
          {-37, 8.533047625744065e-17, 0.9999999999999999,
           8.533047625744065e-17, 8.533047625744063e-17},
          {-30, 9.357622968839299e-14, 0.9999999999999064,
           9.357622968838423e-14, 9.357622968836672e-14},
          {-20, 2.0611536181902037e-09, 0.9999999979388464,
           2.061153613941849e-09, 2.061153605445141e-09},
          {-1, 0.2689414213699951, 0.7310585786300049, 0.19661193324148185,
           0.09085774767294841},
          {0, 0.5, 0.5, 0.25, 0},
          {1, 0.7310585786300049, 0.2689414213699951, 0.19661193324148185,
           -0.09085774767294841},
          {20, 0.9999999979388464, 2.0611536181902037e-09,
           2.061153613941849e-09, -2.061153605445141e-09},
          {30, 0.9999999999999064, 9.357622968839299e-14, 9.357622968838423e-14,
           -9.357622968836672e-14},
          {37, 0.9999999999999999, 8.533047625744065e-17, 8.533047625744065e-17,
           -8.533047625744063e-17},
          {40, 1.0, 4.248354255291589e-18, 4.248354255291589e-18,
           -4.248354255291589e-18},
          {100, 1.0, 3.720075976020836e-44, 3.720075976020836e-44,
           -3.720075976020836e-44},
          {700, 1.0, 9.85967654375977e-305, 9.85967654375977e-305,
           -9.85967654375977e-305},
          {800, 1.0, 0, 0, 0}};
}

// Relative check, absolute below the smallest normal double.
inline void expect_rel(double expected, double actual, double x,
                       double tol = 1e-13) {
  EXPECT_NEAR(
      expected, actual,
      tol * std::fmax(std::fabs(expected), std::numeric_limits<double>::min()))
      << "x = " << x;
}

// First and second derivatives of f at every reference point through
// fvar<double>, fvar<fvar<double>> and fvar<var>, against d1(r) and d2(r).
template <typename F, typename D1, typename D2>
void expect_fwd_derivatives(const F& f, const D1& d1, const D2& d2) {
  using stan::math::fvar;
  using stan::math::var;
  for (const auto& r : refs()) {
    expect_rel(d1(r), f(fvar<double>(r.x, 1)).d_, r.x);

    fvar<fvar<double>> ff(fvar<double>(r.x, 1), fvar<double>(1, 0));
    fvar<fvar<double>> y_ff = f(ff);
    expect_rel(d1(r), y_ff.d_.val_, r.x);
    expect_rel(d2(r), y_ff.d_.d_, r.x);

    fvar<var> fv(r.x, 1);
    fvar<var> y_fv = f(fv);
    expect_rel(d1(r), y_fv.d_.val(), r.x);
    y_fv.d_.grad();
    expect_rel(d2(r), fv.val_.adj(), r.x);
    stan::math::recover_memory();
  }
}

}  // namespace inv_logit_tail_refs

#endif
