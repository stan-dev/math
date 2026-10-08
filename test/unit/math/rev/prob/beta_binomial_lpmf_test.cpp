#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>
#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <vector>

namespace beta_binomial_lpmf_rev_test_internal {
struct TestValue {
  int n;
  int N;
  double alpha;
  double beta;
  double value;
  double grad_log_alpha;  // alpha * d/dalpha
  double grad_log_beta;   // beta * d/dbeta
};

// Computed with mpmath at 80 digits: the value from mp.loggamma, the
// partials from mp.digamma, both checked against mp.diff at 130 digits.
// The shapes are in hex so that they are exact. The gradients are given in
// the log-shape parameterization, which is the one a sampler sees when a
// model puts a prior on log(alpha) or log(concentration).
std::vector<TestValue> testValues = {
    // n = 57, N = 117, alpha = beta = exp(32) / 2
    {57, 117, 0x1.1f43fcc4b662cp+45, 0x1.1f43fcc4b662cp+45, -2.6471538352642870,
     -1.4999999999974545, 1.4999999999981384},
    {400, 1000, 0x1.b48eb57e00000p+44, 0x1.977420dc00000p+42,
     -4.1354072540579752e+2, -4.1081081080252486e+2, 4.1081081078769344e+2},
    {0, 117, 0x1.6345785d8a000p+56, 0x1.0a741a4627800p+58,
     -3.3658802476858363e+1, -2.9249999999999996e+1, 2.9249999999999990e+1},
    {117, 117, 0x1.6345785d8a000p+56, 0x1.0a741a4627800p+58,
     -1.6219644025102715e+2, 8.7749999999999936e+1, -8.7749999999999987e+1},
    // one shape large, the other moderate or small
    {3, 1000000, 0x1.e848000000000p+19, 0x1.bc16d674ec800p+59,
     -4.3238292143128877e+1, 2.9999960000050000, -2.9999989999970000},
    {0, 1000000, 0x1.03caccd133500p+59, 0x1.35c28f5c28f5cp+2,
     -2.8094819745857234e+7, -9.9999999999914529e+5, 5.9751970421184601e+1},
    {1000000, 1000000, 0x1.cfde000000000p+19, 0x1.e000000000000p+3,
     -1.0786783324715004e+1, 7.6922233997281191, -1.0786722596873199e+1},
    {5, 20, 0x1.4000000000000p+3, 0x1.9000000000000p+4, -1.8540068216786033,
     -3.4626330704764243e-1, 5.0911144710080701e-1},
};
}  // namespace beta_binomial_lpmf_rev_test_internal

TEST(ProbDistributionsBetaBinomial, log_shape_gradients) {
  using beta_binomial_lpmf_rev_test_internal::TestValue;
  using beta_binomial_lpmf_rev_test_internal::testValues;
  using stan::math::var;

  for (const TestValue& t : testValues) {
    var alpha = t.alpha;
    var beta = t.beta;
    var lp = stan::math::beta_binomial_lpmf(t.n, t.N, alpha, beta);
    lp.grad();
    EXPECT_NEAR(lp.val(), t.value, 1e-13 * std::max(1.0, std::fabs(t.value)))
        << "n = " << t.n << ", N = " << t.N << ", alpha = " << t.alpha
        << ", beta = " << t.beta;
    // Each partial is a difference of two digamma differences. Its rounding
    // error is a few eps of those terms, which in the log-shape
    // parameterization are of size up to N.
    const double scale = 1e-12 * std::max(1.0, 1e-3 * t.N);
    EXPECT_NEAR(t.alpha * alpha.adj(), t.grad_log_alpha,
                scale * std::max(1.0, std::fabs(t.grad_log_alpha)))
        << "n = " << t.n << ", N = " << t.N << ", alpha = " << t.alpha
        << ", beta = " << t.beta;
    EXPECT_NEAR(t.beta * beta.adj(), t.grad_log_beta,
                scale * std::max(1.0, std::fabs(t.grad_log_beta)))
        << "n = " << t.n << ", N = " << t.N << ", alpha = " << t.alpha
        << ", beta = " << t.beta;
    stan::math::recover_memory();
  }
}

TEST(ProbDistributionsBetaBinomial, log_concentration_gradient) {
  // alpha = beta = exp(lc) / 2, n = 57, N = 117. The
  // gradient in lc goes to 0 like 27 / exp(lc). develop returned 0 or
  // noise (-0.14 at lc = 32) from lc = 20 on, which made the log density a
  // plateau that warmup could not leave.
  using stan::math::var;
  const std::vector<double> lcs = {16.0, 18.0, 20.0, 26.0, 30.0, 32.0};
  const std::vector<double> grads = {
      6.0768267180691362e-6,  8.2241757434662373e-7,  1.1130227121763732e-7,
      2.7589080736553729e-10, 5.0531164031234141e-12, 6.8386493965016465e-13};
  for (size_t i = 0; i < lcs.size(); ++i) {
    var lc = lcs[i];
    var s = stan::math::exp(lc) / 2;
    var lp = stan::math::beta_binomial_lpmf(57, 117, s, s);
    lp.grad();
    EXPECT_NEAR(lc.adj(), grads[i], 1e-12) << "lc = " << lcs[i];
    stan::math::recover_memory();
  }
}
