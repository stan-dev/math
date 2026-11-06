#ifndef TEST_UNIT_MATH_OPENCL_REV_BETA_BINOMIAL_LPMF_REFERENCE_HPP
#define TEST_UNIT_MATH_OPENCL_REV_BETA_BINOMIAL_LPMF_REFERENCE_HPP

#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

namespace beta_binomial_opencl_test_internal {
struct TestValue {
  int n;
  int N;
  double alpha;
  double beta;
  double value;
  double grad_log_alpha;  // alpha * d/dalpha
  double grad_log_beta;   // beta * d/dbeta
};

// The differential OpenCL tests compare OpenCL against the CPU and cannot
// see an error that both share. These are absolute references from mpmath
// at 80 digits: the value from mp.loggamma, the gradients from mp.digamma
// differences, both checked at 130 digits (the gradients with mp.diff). The
// first four rows are as in
// test/unit/math/rev/prob/beta_binomial_lpmf_test.cpp. The shapes are in
// hex so that they are exact. The plain differences of lbeta and digamma
// values keep no correct digits for shapes near 1e15.
// The last row is a small shape where the plain differences are correct.
const std::vector<TestValue> test_values = {
    {57, 117, 0x1.1f43fcc4b662cp+45, 0x1.1f43fcc4b662cp+45, -2.6471538352642870,
     -1.4999999999974545, 1.4999999999981384},
    {400, 1000, 0x1.b48eb57e00000p+44, 0x1.977420dc00000p+42,
     -4.1354072540579752e+2, -4.1081081080252486e+2, 4.1081081078769344e+2},
    {0, 117, 0x1.6345785d8a000p+56, 0x1.0a741a4627800p+58,
     -3.3658802476858363e+1, -2.9249999999999996e+1, 2.9249999999999990e+1},
    {117, 117, 0x1.6345785d8a000p+56, 0x1.0a741a4627800p+58,
     -1.6219644025102715e+2, 8.7749999999999936e+1, -8.7749999999999987e+1},
    // one shape below 10
    {5, 117, 0x1.c6bf526340000p+49, 0x1.35c28f5c28f5cp+2,
     -3.4143943462813216e+3, -1.1199999999999266e+2, 1.5906435196801570e+1},
    {30, 117, 0x1.4000000000000p+1, 0x1.6bcc41e900000p+46,
     -8.2341924150797106e+2, 6.9065498632679532, -2.9999999999966625e+1},
    {5, 20, 0x1.4000000000000p+3, 0x1.9000000000000p+4, -1.8540068216786033,
     -3.4626330704764243e-1, 5.0911144710080701e-1},
};

template <typename T_alpha, typename T_beta>
void expect_reference(const TestValue& t, const stan::math::var& lp,
                      const T_alpha& alpha_adj, const T_beta& beta_adj,
                      const std::string& signature) {
  EXPECT_NEAR(lp.val(), t.value, 1e-12 * std::max(1.0, std::fabs(t.value)))
      << signature << ": n = " << t.n << ", N = " << t.N
      << ", alpha = " << t.alpha << ", beta = " << t.beta;
  EXPECT_NEAR(t.alpha * alpha_adj, t.grad_log_alpha,
              1e-11 * std::max(1.0, std::fabs(t.grad_log_alpha)))
      << signature << ": n = " << t.n << ", N = " << t.N
      << ", alpha = " << t.alpha << ", beta = " << t.beta;
  EXPECT_NEAR(t.beta * beta_adj, t.grad_log_beta,
              1e-11 * std::max(1.0, std::fabs(t.grad_log_beta)))
      << signature << ": n = " << t.n << ", N = " << t.N
      << ", alpha = " << t.alpha << ", beta = " << t.beta;
}
}  // namespace beta_binomial_opencl_test_internal

#endif
