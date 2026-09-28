#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>
#include <gtest/gtest.h>
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

TEST_F(AgradRev, ProbDistributionsBinomial_lccdf_deep_tail) {
  using stan::math::var;
  std::vector<int> n{500};
  var theta = 0.2;
  var lp = stan::math::binomial_lccdf(n, 1000, theta);
  lp.grad();
  EXPECT_NEAR(-227.9265046, lp.val(), 1e-9 * 228.0);
  EXPECT_NEAR(1883.30957, theta.adj(), 1e-8 * 1884.0);
}
