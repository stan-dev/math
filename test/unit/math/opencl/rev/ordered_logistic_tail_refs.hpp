#ifndef TEST_UNIT_MATH_OPENCL_REV_ORDERED_LOGISTIC_TAIL_REFS_HPP
#define TEST_UNIT_MATH_OPENCL_REV_ORDERED_LOGISTIC_TAIL_REFS_HPP

#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <vector>

namespace ordered_logistic_tail_refs {

// Category y, location lambda and cut points (c1, c2); the log probability
// and its derivatives with respect to lambda, c1 and c2.  mpmath, 60 digits,
// rounded to double.  For close cut points the lambda derivative is a
// difference of two values near 1/2; d_lambda_abs_tol > 0 gives an absolute
// tolerance for it.
struct ref {
  int y;
  double lambda;
  double c1;
  double c2;
  double val;
  double d_lambda;
  double d_c1;
  double d_c2;
  double d_lambda_abs_tol;
};

// clang-format off
// NOLINTBEGIN(whitespace/line_length)
inline std::vector<ref> refs() {
  return {
      {2, 1e17, 0, 1, -1e17, -1.0, -0.5819767068693265, 1.5819767068693265, 0},
      {2, -1e17, 0, 1, -1e17, 1.0, -1.5819767068693265, 0.5819767068693265, 0},
      {2, 0.3, 0, 1, -1.416216438741067, 0.09374525535650713, -1.0075341900576675, 0.9137889347011603, 0},
      {1, -40, 0, 1, -4.248354255291589e-18, -4.248354255291589e-18, 4.248354255291589e-18, 0, 0},
      {3, 41, 0, 1, -4.248354255291589e-18, 4.248354255291589e-18, 0, -4.248354255291589e-18, 0},
      {2, 0, 0, 9.313225746154785e-10, -22.18070977791825, 2.3283064365386963e-10, -1073741824.0, 1073741824.0, 1e-15}};
}
// NOLINTEND
// clang-format on

// Relative check, absolute below the smallest normal double, or absolute
// with tolerance abs_tol > 0.
inline void expect_rel(double expected, double actual, const char* what,
                       double lambda, double abs_tol = 0) {
  const double tol = abs_tol > 0
                         ? abs_tol
                         : 1e-13
                               * std::fmax(std::fabs(expected),
                                           std::numeric_limits<double>::min());
  EXPECT_NEAR(expected, actual, tol) << what << " at lambda = " << lambda;
}

}  // namespace ordered_logistic_tail_refs

#endif
