#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>
#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <vector>

namespace {

// beta_lccdf and its shape gradients against a high-precision reference
// (mpmath, 50 digits). The complement 1 - I_y(a, b) is below eps at these
// points, so the value and both gradients are the quantities that a
// linear-space complement cannot represent.
void expect_beta_lccdf(double y, double alpha, double beta, double expected,
                       double dalpha, double dbeta, double rtol) {
  using stan::math::var;
  var a = alpha;
  var b = beta;
  var lp = stan::math::beta_lccdf(y, a, b);
  std::vector<var> x{a, b};
  std::vector<double> g;
  lp.grad(x, g);
  EXPECT_TRUE(std::isfinite(lp.val()))
      << "y=" << y << " a=" << alpha << " b=" << beta;
  EXPECT_TRUE(std::isfinite(g[0]) && std::isfinite(g[1]))
      << "y=" << y << " a=" << alpha << " b=" << beta;
  EXPECT_NEAR(expected, lp.val(), rtol * std::fabs(expected));
  EXPECT_NEAR(dalpha, g[0], rtol * std::fabs(dalpha));
  EXPECT_NEAR(dbeta, g[1], rtol * std::fabs(dbeta));
}

// neg_binomial_lccdf and its alpha and beta gradients against references
// (mpmath, 60 digits, two routes). A value close to 0 is the log of a
// probability that rounds to 1, so the value also has an absolute floor.
void expect_neg_binomial_lccdf(int n, double alpha_in, double beta_in,
                               double expected, double dalpha, double dbeta,
                               double rtol) {
  using stan::math::var;
  var alpha = alpha_in;
  var beta = beta_in;
  var lp = stan::math::neg_binomial_lccdf(n, alpha, beta);
  std::vector<var> x{alpha, beta};
  std::vector<double> g;
  lp.grad(x, g);
  EXPECT_NEAR(expected, lp.val(), std::max(rtol * std::fabs(expected), 1e-15))
      << "n=" << n << " beta=" << beta_in;
  EXPECT_NEAR(dalpha, g[0], rtol * std::fabs(dalpha))
      << "n=" << n << " beta=" << beta_in;
  EXPECT_NEAR(dbeta, g[1], rtol * std::fabs(dbeta))
      << "n=" << n << " beta=" << beta_in;
}

// beta_proportion_lccdf and its mu and kappa gradients, as above
void expect_beta_proportion_lccdf(double y, double mu_in, double kappa_in,
                                  double expected, double dmu, double dkappa,
                                  double rtol) {
  using stan::math::var;
  var mu = mu_in;
  var kappa = kappa_in;
  var lp = stan::math::beta_proportion_lccdf(y, mu, kappa);
  std::vector<var> x{mu, kappa};
  std::vector<double> g;
  lp.grad(x, g);
  EXPECT_NEAR(expected, lp.val(), std::max(rtol * std::fabs(expected), 1e-15))
      << "y=" << y;
  EXPECT_NEAR(dmu, g[0], rtol * std::fabs(dmu)) << "y=" << y;
  EXPECT_NEAR(dkappa, g[1], rtol * std::fabs(dkappa)) << "y=" << y;
}

}  // namespace

TEST_F(AgradRev, ProbDistributionsBeta_lccdf_deep_tail) {
  // moderate shapes, y in the upper tail: the true value is -54, which the
  // complement cannot reach in linear space because eps is e^-36.7
  expect_beta_lccdf(0.9, 50, 50, -54.088884708, 0.595231340835, -1.62672723211,
                    1e-10);
  expect_beta_lccdf(0.8, 50, 50, -25.04429524, 0.481321, -0.937243, 1e-6);
  expect_beta_lccdf(0.758, 352, 590, -316.195574221, 0.708857759363,
                    -0.952707627587, 1e-10);
  // still inside the representable range of the old form, unchanged
  expect_beta_lccdf(0.99, 5, 5, -18.2230299025, 0.737273747257, -4.06065905449,
                    1e-10);
}

TEST_F(AgradRev, ProbDistributionsBeta_proportion_lccdf_deep_tail) {
  using stan::math::var;
  var mu = 0.4;
  var kappa = 900;
  var lp = stan::math::beta_proportion_lccdf(0.758, mu, kappa);
  lp.grad();
  EXPECT_TRUE(std::isfinite(lp.val()));
  EXPECT_TRUE(std::isfinite(mu.adj()) && std::isfinite(kappa.adj()));
  EXPECT_NEAR(-264.2052063, lp.val(), 1e-7 * 264.0);
}

TEST_F(AgradRev, ProbDistributionsNegBinomial_lccdf_deep_tail) {
  using stan::math::var;
  {
    std::vector<int> n{50};
    var alpha = 2;
    var beta = 99;
    var lp = stan::math::neg_binomial_lccdf(n, alpha, beta);
    std::vector<var> x{alpha, beta};
    std::vector<double> g;
    lp.grad(x, g);
    EXPECT_NEAR(-230.9222919, lp.val(), 1e-9 * 231.0);
    EXPECT_NEAR(3.528187864, g[0], 1e-9 * 3.53);
    EXPECT_NEAR(-0.5099009516, g[1], 1e-9 * 0.51);
  }
  {
    std::vector<int> n{50};
    var alpha = 2;
    var beta = 999;
    var lp = stan::math::neg_binomial_lccdf(n, alpha, beta);
    std::vector<var> x{alpha, beta};
    std::vector<double> g;
    lp.grad(x, g);
    EXPECT_NEAR(-348.3452568, lp.val(), 1e-9 * 348.0);
    EXPECT_NEAR(3.5370627, g[0], 1e-9 * 3.54);
    EXPECT_NEAR(-0.05099901827, g[1], 1e-9 * 0.051);
  }
}

TEST_F(AgradRev, ProbDistributionsNegBinomial_lccdf_large_beta) {
  // stan-dev/math#2031: at large beta, p = beta / (beta + 1) rounds towards
  // 1. The complement is formed from 1 / (beta + 1); the density factor of
  // the beta partial must be too, or it is 0 for n >= 1 once p rounds to 1
  // (beta = 1e16, 1e18) and loses digits before that (beta = 1e8). The last
  // two rows are controls on both sides of beta = 1.
  const double rtol = 1e-12;
  // the reported case, log1m_exp(-alpha * log1p(1 / beta))
  expect_neg_binomial_lccdf(0, 1.0, 1e18, -41.446531673892821, 1.0,
                            -1.0000000000000001e-18, rtol);
  expect_neg_binomial_lccdf(1, 2.0, 1e18, -81.794451059117534,
                            0.83333333333333337, -2.0000000000000001e-18, rtol);
  expect_neg_binomial_lccdf(1, 2.0, 1e16, -72.584110687141347,
                            0.83333333333333326, -1.9999999999999997e-16, rtol);
  expect_neg_binomial_lccdf(3, 2.0, 1e8, -72.073285111375355,
                            1.2833333253333334, -3.9999999520000002e-08, rtol);
  expect_neg_binomial_lccdf(5, 0.5, 1e-3, -0.08928900804336741,
                            0.48141786839236278, -46.496417635081194, rtol);
  expect_neg_binomial_lccdf(5, 3.0, 2.0, -3.9290859049832054,
                            0.8864463205559171, -1.7364341085271318, rtol);
}

TEST_F(AgradRev, ProbDistributionsBeta_lccdf_small_y) {
  // The complement 1 - I_y(a, b) and its shape gradients must be evaluated
  // from y, not from 1 - y, which is close to 1 for small y: from 1 - y the
  // alpha gradient at y = 1e-10 had a relative error of 4e-8.
  expect_beta_lccdf(1e-10, 0.5, 3, -1.8750175782197275e-05,
                    0.00041174242508107103, -3.3820441414792342e-06, 1e-11);
  expect_beta_lccdf(1e-6, 0.5, 0.5, -0.00063682260695097634,
                    0.0091917773883985632, -0.00088310453739762977, 1e-11);
  // y <= 1/2 with the complement far below eps
  expect_beta_lccdf(0.3, 1, 200, -71.334988787746468, 4.6855356557158476,
                    -0.35667494393873234, 1e-12);
  expect_beta_lccdf(0.45, 50, 450, -150.59662983495602, 1.5162204045412182,
                    -0.49492426185967869, 1e-12);
}

TEST_F(AgradRev, ProbDistributionsBeta_proportion_lccdf_small_y) {
  // as for beta_lccdf; the value -8.4e-23 is the log of a probability that
  // rounds to 1
  expect_beta_proportion_lccdf(1e-8, 0.3, 10, -8.3999996220000247e-23,
                               1.4955371148778298e-20, 4.1682780319251551e-22,
                               1e-11);
  expect_beta_proportion_lccdf(0.4, 0.02, 500, -215.48431197855072,
                               1770.6978109683873, -0.42187260226452955, 1e-12);
}

TEST_F(AgradRev, ProbDistributionsNegBinomial_lccdf_small_beta) {
  // beta < 1, so p = beta / (beta + 1) < 1/2: the complement and the alpha
  // gradient must be evaluated from p, not from 1 / (beta + 1), which is
  // close to 1: from 1 / (beta + 1) the alpha gradient at beta = 1e-9 had a
  // relative error of 4e-8. The last two rows have the complement far below
  // eps.
  expect_neg_binomial_lccdf(5, 2.0, 1e-6, -2.0999888000598496e-11,
                            2.6717432945696996e-10, -4.1999664002393978e-05,
                            1e-11);
  expect_neg_binomial_lccdf(5, 0.5, 1e-9, -8.56075085052458e-05,
                            0.0016237738030832265, -42805.586280793519, 1e-11);
  expect_neg_binomial_lccdf(100, 1.0, 0.5, -40.951975918924603,
                            4.1179071920897643, -67.333333333333329, 1e-12);
  expect_neg_binomial_lccdf(300, 2.0, 0.2, -50.943700316621204,
                            3.5140890478701055, -246.74809989142236, 1e-12);
}

TEST_F(AgradRev, ProbDistributionsBinomial_lccdf_deep_tail) {
  using stan::math::var;
  std::vector<int> n{500};
  var theta = 0.2;
  var lp = stan::math::binomial_lccdf(n, 1000, theta);
  lp.grad();
  EXPECT_NEAR(-227.9265046, lp.val(), 1e-9 * 228.0);
  EXPECT_NEAR(1883.30957, theta.adj(), 1e-8 * 1884.0);
}
