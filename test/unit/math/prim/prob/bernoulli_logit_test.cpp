#include <stan/math/prim.hpp>
#include <stan/math/rev.hpp>
#include <test/unit/math/prim/prob/vector_rng_test_helper.hpp>
#include <test/unit/math/prim/prob/VectorIntRNGTestRig.hpp>
#include <test/unit/math/prim/prob/util.hpp>
#include <boost/math/distributions.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <vector>

class BernoulliLogitTestRig : public VectorIntRNGTestRig {
 public:
  BernoulliLogitTestRig()
      : VectorIntRNGTestRig(10000, 10, {0, 1},
                            {-5.7, -1.0, 0.0, 0.2, 1.0, 10.0}, {-3, -2, 0, 1},
                            {}, {}) {}

  template <typename T1, typename T2, typename T3, typename T_rng>
  auto generate_samples(const T1& t, const T2&, const T3&, T_rng& rng) const {
    return stan::math::bernoulli_logit_rng(t, rng);
  }

  template <typename T1>
  double pmf(int y, T1 t, double, double) const {
    return std::exp(stan::math::bernoulli_logit_lpmf(y, t));
  }
};

TEST(ProbDistributionsBernoulliLogit, errorCheck) {
  check_dist_throws_all_types(BernoulliLogitTestRig());
}

TEST(ProbDistributionsBernoulliLogit, distributionCheck) {
  check_counts_real(BernoulliLogitTestRig());
}

TEST(ProbDistributionsBernoulliLogit, error_check) {
  boost::random::mt19937 rng;

  EXPECT_NO_THROW(stan::math::bernoulli_logit_rng(-3.5, rng));
  EXPECT_THROW(
      stan::math::bernoulli_logit_rng(stan::math::positive_infinity(), rng),
      std::domain_error);
}

TEST(ProbDistributionsBernoulliLogit, logitChiSquareGoodnessFitTest) {
  boost::random::mt19937 rng;
  // number of samples
  int N = 10000;

  // logit-transformed probability
  double parameter = -0.5;
  // actual probability
  double prob = stan::math::inv_logit(-0.5);

  std::vector<double> expected;
  expected.push_back(N * (1 - prob));
  expected.push_back(N * prob);

  std::vector<int> counts(2);
  for (int i = 0; i < N; ++i) {
    ++counts[stan::math::bernoulli_logit_rng(parameter, rng)];
  }

  assert_chi_squared(counts, expected, 1e-6);
}

TEST(ProbDistributionsBernoulliLogit, cutoff) {
  double cutoff = 20;
  for (int n = 0; n <= 1; ++n) {
    for (int sign : {-1, 1}) {
      double before_cutoff
          = stan::math::bernoulli_logit_lpmf(n, sign * cutoff - 1e-14);
      double after_cutoff
          = stan::math::bernoulli_logit_lpmf(n, sign * cutoff + 1e-14);
      double relative_error_at_cutoff = log(before_cutoff / after_cutoff);
      EXPECT_NEAR(relative_error_at_cutoff, 0, 1e-8)
          << "bernoulli_logit_lpmf changes too much around cutoff for n = " << n
          << ", cutoff = " << (sign * cutoff)
          << ", value at cutoff - 1e-14: " << before_cutoff
          << ", value at cutoff + 1e-14: " << after_cutoff;
    }
  }
}

TEST(ProbDistributionsBernoulliLogit, cutoff_partials_sign) {
  // For ntheta > cutoff the value is -exp(-ntheta), so the derivative with
  // respect to theta is signs * exp(-ntheta): positive for y = 1, negative
  // for y = 0. Regression test: the (ntheta > cutoff) partials branch used
  // to return -exp_m_ntheta without the signs factor, which has the wrong
  // sign whenever y = 1 and theta > cutoff (and, symmetrically, for y = 0
  // with theta < -cutoff).
  using stan::math::var;
  const double cutoff = 20;
  const double h = 1e-3;  // keeps both FD points inside the same branch
  for (int n = 0; n <= 1; ++n) {
    const double theta = (n == 1 ? 1.0 : -1.0) * (cutoff + 5);
    const double signs = 2 * n - 1;
    const double expected = signs * std::exp(-signs * theta);

    // finite-difference reference from the double implementation
    const double fd = (stan::math::bernoulli_logit_lpmf(n, theta + h)
                       - stan::math::bernoulli_logit_lpmf(n, theta - h))
                      / (2 * h);
    EXPECT_NEAR(expected, fd, 1e-6 * std::fabs(fd))
        << "analytic derivative disagrees with finite differences for n = " << n;

    // autodiff gradient
    var theta_v = theta;
    var lp = stan::math::bernoulli_logit_lpmf(n, theta_v);
    lp.grad();
    EXPECT_NEAR(theta_v.adj(), expected, 1e-9 * std::fabs(expected))
        << "gradient has the wrong sign or magnitude for n = " << n
        << ", theta = " << theta << ", adj = " << theta_v.adj()
        << ", expected = " << expected;
    stan::math::recover_memory();
  }
}
