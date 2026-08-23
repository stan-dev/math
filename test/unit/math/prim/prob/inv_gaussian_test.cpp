#include <stan/math/prim.hpp>
#include <test/unit/math/prim/prob/vector_rng_test_helper.hpp>
#include <test/unit/math/prim/prob/util.hpp>
#include <gtest/gtest.h>
#include <boost/math/distributions.hpp>
#include <cmath>
#include <limits>
#include <vector>

class InvGaussianTestRig : public VectorRealRNGTestRig {
 public:
  InvGaussianTestRig()
      : VectorRealRNGTestRig(10000, 10, {0.1, 0.5, 1.0, 2.5, 4.0}, {1, 2, 3, 4},
                             {0.0, -0.1, -1.0}, {0, -1, -2},
                             {0.2, 1.0, 3.0, 7.5}, {1, 2, 3, 8},
                             {0.0, -0.1, -1.0}, {0, -1, -2}) {}

  template <typename T1, typename T2, typename T3, typename T_rng>
  auto generate_samples(const T1& mu, const T2& lambda, const T3& unused,
                        T_rng& rng) const {
    return stan::math::inv_gaussian_rng(mu, lambda, rng);
  }

  std::vector<double> generate_quantiles(double mu, double lambda,
                                         double unused) const {
    std::vector<double> quantiles;
    double K = stan::math::round(2 * std::pow(N_, 0.4));
    boost::math::inverse_gaussian_distribution<> dist(mu, lambda);
    for (int i = 1; i < K; ++i) {
      double frac = i / K;
      quantiles.push_back(quantile(dist, frac));
    }
    quantiles.push_back(std::numeric_limits<double>::max());
    return quantiles;
  }
};

TEST(ProbDistributionsInvGaussian, errorCheck) {
  check_dist_throws_all_types(InvGaussianTestRig());
}

TEST(ProbDistributionsInvGaussian, distributionCheck) {
  check_quantiles_real_real(InvGaussianTestRig());
}

TEST(ProbDistributionsInvGaussian, error_check) {
  boost::random::mt19937 rng;
  EXPECT_NO_THROW(stan::math::inv_gaussian_rng(1.0, 2.0, rng));
  EXPECT_THROW(stan::math::inv_gaussian_rng(0.0, 2.0, rng), std::domain_error);
  EXPECT_THROW(stan::math::inv_gaussian_rng(-1.0, 2.0, rng), std::domain_error);
  EXPECT_THROW(stan::math::inv_gaussian_rng(1.0, 0.0, rng), std::domain_error);
  EXPECT_THROW(stan::math::inv_gaussian_rng(1.0, -1.0, rng), std::domain_error);
  EXPECT_THROW(
      stan::math::inv_gaussian_rng(stan::math::positive_infinity(), 2.0, rng),
      std::domain_error);
  EXPECT_THROW(
      stan::math::inv_gaussian_rng(1.0, stan::math::positive_infinity(), rng),
      std::domain_error);
}

TEST(ProbDistributionsInvGaussian, rngStableForLargeMuOverLambda) {
  boost::random::mt19937 rng(1234);
  for (double mu : {1.0, 1e3, 1e6, 1e9, 1e12}) {
    for (double lambda : {1e-8, 1e-4, 1.0, 1e4}) {
      for (int i = 0; i < 2000; ++i) {
        double d = stan::math::inv_gaussian_rng(mu, lambda, rng);
        ASSERT_TRUE(std::isfinite(d)) << "mu=" << mu << " lambda=" << lambda;
        ASSERT_GT(d, 0.0) << "mu=" << mu << " lambda=" << lambda;
      }
    }
  }
}

// The sample variance of an inverse Gaussian is heavy tailed, so its
// tolerance is wider than the mean's.
TEST(ProbDistributionsInvGaussian, rngMomentsAtLargeMu) {
  boost::random::mt19937 rng(4321);
  const double mu = 1e3;
  const double lambda = 1e3;
  const int N = 200000;
  double sum = 0;
  double sum_sq = 0;
  for (int i = 0; i < N; ++i) {
    double d = stan::math::inv_gaussian_rng(mu, lambda, rng);
    sum += d;
    sum_sq += d * d;
  }
  double mean = sum / N;
  double var = sum_sq / N - mean * mean;
  EXPECT_NEAR(mu, mean, 0.05 * mu);
  EXPECT_NEAR(mu * mu * mu / lambda, var, 0.15 * mu * mu * mu / lambda);
}

// reference values from mpmath at 60 digits

TEST(ProbDistributionsInvGaussian, values) {
  using stan::math::inv_gaussian_cdf;
  using stan::math::inv_gaussian_lccdf;
  using stan::math::inv_gaussian_lcdf;
  using stan::math::inv_gaussian_lpdf;

  EXPECT_FLOAT_EQ(-2.479180611448965157415, inv_gaussian_lpdf(1.2, 0.5, 2.0));
  EXPECT_FLOAT_EQ(-2.391593703832052124008, inv_gaussian_lpdf(0.3, 1.0, 5.0));
  EXPECT_FLOAT_EQ(0.9815922531042920910742, inv_gaussian_cdf(1.2, 0.5, 2.0));
  EXPECT_FLOAT_EQ(0.003359190912064955560317, inv_gaussian_cdf(0.3, 1.0, 5.0));
  EXPECT_FLOAT_EQ(-0.01857927772712011847827, inv_gaussian_lcdf(1.2, 0.5, 2.0));
  EXPECT_FLOAT_EQ(-5.696055133984662592139, inv_gaussian_lcdf(0.3, 1.0, 5.0));
  EXPECT_FLOAT_EQ(-3.9949836760335225255, inv_gaussian_lccdf(1.2, 0.5, 2.0));
  EXPECT_FLOAT_EQ(-0.003364845660995599504277,
                  inv_gaussian_lccdf(0.3, 1.0, 5.0));
}

TEST(ProbDistributionsInvGaussian, boundaries) {
  using stan::math::inv_gaussian_cdf;
  using stan::math::inv_gaussian_lccdf;
  using stan::math::inv_gaussian_lcdf;
  using stan::math::inv_gaussian_lpdf;
  double inf = std::numeric_limits<double>::infinity();

  EXPECT_FLOAT_EQ(-inf, inv_gaussian_lpdf(0.0, 1.0, 2.0));
  EXPECT_FLOAT_EQ(-inf, inv_gaussian_lcdf(0.0, 1.0, 2.0));
  EXPECT_FLOAT_EQ(0.0, inv_gaussian_lccdf(0.0, 1.0, 2.0));
  EXPECT_FLOAT_EQ(0.0, inv_gaussian_cdf(0.0, 1.0, 2.0));
  EXPECT_FLOAT_EQ(0.0, inv_gaussian_lcdf(inf, 1.0, 2.0));
  EXPECT_FLOAT_EQ(-inf, inv_gaussian_lccdf(inf, 1.0, 2.0));
  EXPECT_FLOAT_EQ(1.0, inv_gaussian_cdf(inf, 1.0, 2.0));
}

// 2 lambda / mu is past the overflow point of exp for all of these.
TEST(ProbDistributionsInvGaussian, largeExpFactor) {
  using stan::math::inv_gaussian_lccdf;
  using stan::math::inv_gaussian_lcdf;

  // 2 lambda / mu = 2000
  EXPECT_FLOAT_EQ(-0.68061354470344461635, inv_gaussian_lcdf(0.1, 0.1, 100.0));
  EXPECT_FLOAT_EQ(-0.70583990450544626626, inv_gaussian_lccdf(0.1, 0.1, 100.0));
  EXPECT_FLOAT_EQ(-671.8777811569157, inv_gaussian_lccdf(0.3, 0.1, 100.0));
  EXPECT_FLOAT_EQ(-4057.1238032952494, inv_gaussian_lccdf(1.0, 0.1, 100.0));
  // 2 lambda / mu = 200
  EXPECT_FLOAT_EQ(-805.707747545284162, inv_gaussian_lccdf(0.5, 0.1, 50.0));
}

// erfc underflows to zero here; the asymptotic log_Phi branch carries these.
TEST(ProbDistributionsInvGaussian, deepLowerTail) {
  using stan::math::inv_gaussian_lcdf;

  EXPECT_FLOAT_EQ(-906.524245083721553,
                  inv_gaussian_lcdf(0.05, 1.0, 100.0));  // z1 = -42.5
  EXPECT_FLOAT_EQ(-849.06768101617034,
                  inv_gaussian_lcdf(0.1, 0.75, 225.0));  // z1 = -41.1
  EXPECT_FLOAT_EQ(-5333.88923087778086,
                  inv_gaussian_lcdf(0.02, 0.75, 225.0));  // z1 = -103.2
  EXPECT_FLOAT_EQ(-4905.33096155861673,
                  inv_gaussian_lcdf(0.01, 1.0, 100.0));  // z1 = -99.0
}

// Once Phi(z1) rounds to exactly one the log CDF lands a few 1e-20 above
// zero. That is bounded by the rounding of Phi and vanishes under exp, so it
// is asserted at that scale.
TEST(ProbDistributionsInvGaussian, probabilityNeverExceedsOne) {
  using stan::math::inv_gaussian_lccdf;
  using stan::math::inv_gaussian_lcdf;
  for (double mu : {1e-3, 1e-2, 0.1, 1.0, 10.0, 1e3}) {
    for (double lambda : {1.0, 1e3, 1e8, 1e14, 1e16}) {
      for (double rel : {1e-3, 0.1, 0.5, 1.0, 2.0, 10.0, 1e3}) {
        double y = mu * rel;
        double lf = inv_gaussian_lcdf(y, mu, lambda);
        double ls = inv_gaussian_lccdf(y, mu, lambda);
        EXPECT_LE(lf, 1e-15);
        EXPECT_LE(ls, 1e-15);
        EXPECT_LE(std::exp(lf), 1.0);
        EXPECT_LE(std::exp(ls), 1.0);
      }
    }
  }
  // against mpmath at 60 digits
  EXPECT_NEAR(-0.693147179298379049, inv_gaussian_lcdf(0.1, 0.1, 1e16), 1e-12);
  EXPECT_FLOAT_EQ(-4.04999999999999999e17, inv_gaussian_lcdf(1e-4, 1e-3, 1e14));
  EXPECT_FLOAT_EQ(-4.99000499999999997e17, inv_gaussian_lcdf(1e-4, 1e-1, 1e14));
}

// At y == mu the CDF is a fixed quantity independent of the scale of mu, so
// all of these must agree.
TEST(ProbDistributionsInvGaussian, medianIsScaleInvariant) {
  using stan::math::inv_gaussian_lccdf;
  using stan::math::inv_gaussian_lcdf;
  for (double lambda_over_mu : {1e2, 1e6, 1e11, 1e13}) {
    double ref_lcdf = 0;
    double ref_lccdf = 0;
    bool first = true;
    for (double mu : {1e-3, 1e-2, 0.1, 1.0, 10.0}) {
      double lambda = lambda_over_mu * mu;
      double a = inv_gaussian_lcdf(mu, mu, lambda);
      double b = inv_gaussian_lccdf(mu, mu, lambda);
      if (first) {
        ref_lcdf = a;
        ref_lccdf = b;
        first = false;
      } else {
        EXPECT_NEAR(ref_lcdf, a, 1e-13);
        EXPECT_NEAR(ref_lccdf, b, 1e-13);
      }
    }
  }
  // F(mu) for a very large shape approaches 1/2 from above
  EXPECT_NEAR(-0.693147181822, inv_gaussian_lccdf(1e-3, 1e-3, 1e14), 1e-11);
  EXPECT_NEAR(-0.693147180686, inv_gaussian_lccdf(1e-3, 1e-3, 1e16), 1e-11);
}

TEST(ProbDistributionsInvGaussian, cdfCcdfSumToOne) {
  using stan::math::inv_gaussian_lccdf;
  using stan::math::inv_gaussian_lcdf;
  std::vector<double> mus{0.5, 1.0, 0.3, 2.0};
  std::vector<double> lambdas{2.0, 5.0, 10.0, 0.7};
  std::vector<double> ys{0.1, 0.3, 0.8, 1.5, 4.0};
  for (double mu : mus) {
    for (double lambda : lambdas) {
      for (double y : ys) {
        double f = std::exp(inv_gaussian_lcdf(y, mu, lambda));
        double s = std::exp(inv_gaussian_lccdf(y, mu, lambda));
        EXPECT_NEAR(1.0, f + s, 1e-12);
      }
    }
  }
}

TEST(ProbDistributionsInvGaussian, vectorMatchesScalarSum) {
  using stan::math::inv_gaussian_lccdf;
  using stan::math::inv_gaussian_lcdf;
  using stan::math::inv_gaussian_lpdf;
  std::vector<double> y{0.2, 0.5, 1.0, 2.0, 5.0};
  double sum_lpdf = 0;
  double sum_lcdf = 0;
  double sum_lccdf = 0;
  for (double yi : y) {
    sum_lpdf += inv_gaussian_lpdf(yi, 1.3, 3.0);
    sum_lcdf += inv_gaussian_lcdf(yi, 1.3, 3.0);
    sum_lccdf += inv_gaussian_lccdf(yi, 1.3, 3.0);
  }
  EXPECT_FLOAT_EQ(sum_lpdf, inv_gaussian_lpdf(y, 1.3, 3.0));
  EXPECT_FLOAT_EQ(sum_lcdf, inv_gaussian_lcdf(y, 1.3, 3.0));
  EXPECT_FLOAT_EQ(sum_lccdf, inv_gaussian_lccdf(y, 1.3, 3.0));
}

// check helper functions; tolerances scale with each value's magnitude
TEST(ProbDistributionsInvGaussian, internalLogPhi) {
  using stan::math::internal::log_Phi;

  // erfc branch, up to the switch at z = -30
  EXPECT_NEAR(-0.6931471805599453094, log_Phi(0.0), 1e-15);
  EXPECT_NEAR(-1.841021645009263506, log_Phi(-1.0), 1e-14);
  EXPECT_NEAR(-15.06499839398872574, log_Phi(-5.0), 1e-13);
  EXPECT_NEAR(-53.23128515051247058, log_Phi(-10.0), 1e-12);
  EXPECT_NEAR(-203.9171553710972639, log_Phi(-20.0), 1e-11);
  EXPECT_NEAR(-451.3229124585286345, log_Phi(-29.9), 1e-11);
  // asymptotic branch
  EXPECT_NEAR(-454.3212439563431971, log_Phi(-30.0), 1e-11);
  EXPECT_NEAR(-457.3295644163822579, log_Phi(-30.1), 1e-11);
  // beyond the point where erfc underflows to zero (z < -37.5)
  EXPECT_NEAR(-745.6952702904110813, log_Phi(-38.5), 1e-10);
  EXPECT_NEAR(-804.6084420137537882, log_Phi(-40.0), 1e-10);
  EXPECT_NEAR(-1805.013560680567139, log_Phi(-60.0), 1e-9);
  EXPECT_NEAR(-11255.92961826680818, log_Phi(-150.0), 1e-8);
  EXPECT_NEAR(-500007.8266948121843, log_Phi(-1000.0), 1e-6);
  // upper tail saturates at log(1) = 0
  EXPECT_NEAR(-2.866516129637635934e-7, log_Phi(5.0), 1e-16);
  EXPECT_FLOAT_EQ(0.0, log_Phi(50.0));

  double inf = std::numeric_limits<double>::infinity();
  EXPECT_FLOAT_EQ(0.0, log_Phi(inf));
  EXPECT_FLOAT_EQ(-inf, log_Phi(-inf));
  EXPECT_TRUE(std::isnan(log_Phi(std::numeric_limits<double>::quiet_NaN())));
}

// The true slope at z = -30 is about 30, so across a 2e-10 interval the
// honest change is about 6e-9.
TEST(ProbDistributionsInvGaussian, internalLogPhiBranchContinuity) {
  using stan::math::internal::log_Phi;
  double eps = 1e-10;
  double step = log_Phi(-30.0 + eps) - log_Phi(-30.0 - eps);
  EXPECT_LT(std::fabs(step), 1e-7);
}

TEST(ProbDistributionsInvGaussian, errors) {
  using stan::math::inv_gaussian_lpdf;
  double inf = std::numeric_limits<double>::infinity();

  EXPECT_THROW(inv_gaussian_lpdf(-1.0, 1.0, 2.0), std::domain_error);
  EXPECT_THROW(inv_gaussian_lpdf(1.0, 0.0, 2.0), std::domain_error);
  EXPECT_THROW(inv_gaussian_lpdf(1.0, -1.0, 2.0), std::domain_error);
  EXPECT_THROW(inv_gaussian_lpdf(1.0, inf, 2.0), std::domain_error);
  EXPECT_THROW(inv_gaussian_lpdf(1.0, 1.0, 0.0), std::domain_error);
  EXPECT_THROW(inv_gaussian_lpdf(1.0, 1.0, -1.0), std::domain_error);
  EXPECT_THROW(inv_gaussian_lpdf(1.0, 1.0, inf), std::domain_error);

  std::vector<double> y{1.0, 2.0};
  std::vector<double> mu{1.0, 2.0, 3.0};
  EXPECT_THROW(inv_gaussian_lpdf(y, mu, 1.0), std::invalid_argument);
}
