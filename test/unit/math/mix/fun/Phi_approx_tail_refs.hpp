#ifndef TEST_UNIT_MATH_MIX_FUN_PHI_APPROX_TAIL_REFS_HPP
#define TEST_UNIT_MATH_MIX_FUN_PHI_APPROX_TAIL_REFS_HPP

#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <vector>

// First and second derivatives of Phi_approx(x) = inv_logit(c3 x^3 + c1 x),
// computed with mpmath at 60 digits from the double values of
// PHI_APPROX_CUBIC_COEF (c3) and PHI_APPROX_LINEAR_COEF (c1). At x = 40 and
// -40 the true values, about 3e-1987 and 1e-1984, underflow to 0.
namespace Phi_approx_tail_refs {

inline std::vector<double> args() {
  return {-40, -10, -8, -7, -1, 0, 1, 7, 8, 10, 40};
}

inline std::vector<double> derivs() {
  return {0.0,
          5.934767576595741e-37,
          8.691773531236104e-21,
          5.126988948527858e-15,
          0.2415281648469019,
          0.399400718560196,
          0.2415281648469019,
          5.126988948527858e-15,
          8.691773531236104e-21,
          5.934767576595741e-37,
          0.0};
}

inline std::vector<double> second_derivs() {
  return {0.0,
          1.3401219329217354e-35,
          1.297010750153944e-19,
          6.010339719762636e-14,
          0.24180176878701165,
          0.0,
          -0.24180176878701165,
          -6.010339719762636e-14,
          -1.297010750153944e-19,
          -1.3401219329217354e-35,
          0.0};
}

// Relative check, absolute below the smallest normal double.
inline void expect_rel(double expected, double actual, double x,
                       double tol = 1e-13) {
  EXPECT_NEAR(
      expected, actual,
      tol * std::fmax(std::fabs(expected), std::numeric_limits<double>::min()))
      << "x = " << x;
}

}  // namespace Phi_approx_tail_refs

#endif
