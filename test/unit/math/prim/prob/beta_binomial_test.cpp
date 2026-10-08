#include <stan/math/prim.hpp>
#include <test/unit/math/prim/prob/vector_rng_test_helper.hpp>
#include <test/unit/math/prim/prob/VectorIntRNGTestRig.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/math/distributions.hpp>
#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

class BetaBinomialTestRig : public VectorIntRNGTestRig {
 public:
  BetaBinomialTestRig()
      : VectorIntRNGTestRig(10000, 10, {0, 1, 2, 3, 4, 5, 6}, {}, {0, 1, 3, 8},
                            {}, {-1, -5, -7}, {0.1, 1.7, 3.99}, {1, 2, 3},
                            {-2.1, -0.5, 0.0}, {-3, -1, 0}, {0.1, 1.1, 4.99},
                            {1, 2, 3}, {-3.0, -2.0, 0.0}, {-3, -1, 0}) {}

  template <typename T1, typename T2, typename T3, typename T_rng>
  auto generate_samples(const T1& N, const T2& alpha, const T3& beta,
                        T_rng& rng) const {
    return stan::math::beta_binomial_rng(N, alpha, beta, rng);
  }

  template <typename T1>
  double pmf(int y, T1 N, double alpha, double beta) const {
    if (y <= N) {
      return std::exp(stan::math::beta_binomial_lpmf(y, N, alpha, beta));
    } else {
      return 0.0;
    }
  }
};

TEST(ProbDistributionsBetaBinomial, errorCheck) {
  check_dist_throws_int_first_argument(BetaBinomialTestRig());
}

TEST(ProbDistributionsBetaBinomial, distributionCheck) {
  check_counts_int_real_real(BetaBinomialTestRig());
}

TEST(ProbDistributionBetaBinomial, error_check) {
  boost::random::mt19937 rng;
  EXPECT_NO_THROW(stan::math::beta_binomial_rng(4, 0.6, 2.0, rng));

  EXPECT_THROW(stan::math::beta_binomial_rng(-4, 0.6, 2, rng),
               std::domain_error);
  EXPECT_THROW(stan::math::beta_binomial_rng(4, -0.6, 2, rng),
               std::domain_error);
  EXPECT_THROW(stan::math::beta_binomial_rng(4, 0.6, -2, rng),
               std::domain_error);
  EXPECT_THROW(
      stan::math::beta_binomial_rng(4, stan::math::positive_infinity(), 2, rng),
      std::domain_error);
  EXPECT_THROW(stan::math::beta_binomial_rng(
                   4, 0.6, stan::math::positive_infinity(), rng),
               std::domain_error);
}

namespace beta_binomial_test_internal {
struct TestValue {
  int n;
  int N;
  double alpha;
  double beta;
  double value;
};

// Log pmf computed with mpmath at 80 digits as lchoose(N, n)
// + lbeta(n + alpha, N - n + beta) - lbeta(alpha, beta) with mp.loggamma,
// and checked at 130 digits. The shapes are written in hex so that they are
// exact. The first three are n = 57, N = 117, alpha = beta = exp(lc) / 2
// for lc = 32, 36, 42, where the plain lbeta difference was wrong by
// 4.5e-3, 9.8e-2 and 81. The next three are n = 500, N = 1000,
// alpha = beta = 0.5e10, 0.5e13, 0.5e19, where the plain lbeta difference
// was wrong by 5.3e-8, 2.8e-4 and 693. In the last two points one or both
// shapes are small, where the lbeta form is the accurate one.
std::vector<TestValue> testValues = {
    {57, 117, 0x1.1f43fcc4b662cp+45, 0x1.1f43fcc4b662cp+45,
     -2.6471538352642870},
    {57, 117, 0x1.ea215a1d20d76p+50, 0x1.ea215a1d20d76p+50,
     -2.6471538352636157},
    {57, 117, 0x1.8232558201159p+59, 0x1.8232558201159p+59,
     -2.6471538352636032},
    {500, 1000, 0x1.2a05f20000000p+32, 0x1.2a05f20000000p+32,
     -3.6799190420941268},
    {500, 1000, 0x1.2309ce5400000p+42, 0x1.2309ce5400000p+42,
     -3.6799189921441293},
    {500, 1000, 0x1.158e460913d00p+62, 0x1.158e460913d00p+62,
     -3.6799189920941293},
    {400, 1000, 0x1.b48eb57e00000p+44, 0x1.977420dc00000p+42,
     -4.1354072540579752e+2},
    {0, 117, 0x1.6345785d8a000p+56, 0x1.0a741a4627800p+58,
     -3.3658802476858363e+1},
    {117, 117, 0x1.6345785d8a000p+56, 0x1.0a741a4627800p+58,
     -1.6219644025102715e+2},
    {3, 1000000, 0x1.e848000000000p+19, 0x1.bc16d674ec800p+59,
     -4.3238292143128877e+1},
    {1000000, 1000000, 0x1.cfde000000000p+19, 0x1.e000000000000p+3,
     -1.0786783324715004e+1},
    {5, 20, 0x1.4000000000000p+3, 0x1.9000000000000p+4, -1.8540068216786033},
};
}  // namespace beta_binomial_test_internal

TEST(ProbDistributionsBetaBinomial, large_shapes) {
  using beta_binomial_test_internal::TestValue;
  using beta_binomial_test_internal::testValues;
  using stan::math::beta_binomial_lpmf;

  for (const TestValue& t : testValues) {
    const double tol = 1e-13 * std::max(1.0, std::fabs(t.value));
    EXPECT_NEAR(beta_binomial_lpmf(t.n, t.N, t.alpha, t.beta), t.value, tol)
        << "n = " << t.n << ", N = " << t.N << ", alpha = " << t.alpha
        << ", beta = " << t.beta;
  }
}

TEST(ProbDistributionsBetaBinomial, binomial_limit) {
  // The difference to the binomial limit tends to 0
  // from below as the shapes grow; at lc = 42 it is -3.1e-17
  using stan::math::beta_binomial_lpmf;
  using stan::math::binomial_lpmf;
  const double s = 0x1.8232558201159p+59;  // exp(42) / 2
  const double diff
      = beta_binomial_lpmf(57, 117, s, s) - binomial_lpmf(57, 117, 0.5);
  EXPECT_LE(diff, 1e-13);
  EXPECT_GE(diff, -1e-13);
}
