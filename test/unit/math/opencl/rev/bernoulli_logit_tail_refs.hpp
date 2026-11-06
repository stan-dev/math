#ifndef TEST_UNIT_MATH_OPENCL_REV_BERNOULLI_LOGIT_TAIL_REFS_HPP
#define TEST_UNIT_MATH_OPENCL_REV_BERNOULLI_LOGIT_TAIL_REFS_HPP

#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <vector>

namespace bernoulli_logit_tail_refs {

// Outcome y and logit theta; the log probability log(inv_logit(t)), with
// t = (2y - 1) theta, and its derivative (2y - 1) inv_logit(-t) with
// respect to theta.  mpmath, 60 digits, rounded to double.
struct ref {
  int y;
  double theta;
  double val;
  double d_theta;
};

inline std::vector<ref> refs() {
  return {{1, 21.0, -7.582560425037146e-10, 7.582560422162384e-10},
          {0, 21.0, -21.000000000758256, -0.999999999241744},
          {1, -25.0, -25.000000000013888, 0.9999999999861121},
          {0, -25.0, -1.3887943864867583e-11, -1.3887943864771146e-11},
          {1, 30.0, -9.357622968839737e-14, 9.357622968839299e-14},
          {0, 30.0, -30.000000000000092, -0.9999999999999064},
          {1, 0.3, -0.5543552444685271, 0.425557483188341}};
}

// Relative check, absolute below the smallest normal double.
inline void expect_rel(double expected, double actual, const char* what,
                       const ref& r) {
  EXPECT_NEAR(
      expected, actual,
      1e-13
          * std::fmax(std::fabs(expected), std::numeric_limits<double>::min()))
      << what << " at y = " << r.y << ", theta = " << r.theta;
}

}  // namespace bernoulli_logit_tail_refs

#endif
