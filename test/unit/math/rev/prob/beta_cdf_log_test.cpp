#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <vector>

namespace {

// Shape gradients of beta_lcdf against a high-precision reference: the
// reference is (d/da I) / I and (d/db I) / I from mpmath at 80 digits.
void expect_beta_lcdf_grads(double y, double alpha, double beta, double dalpha,
                            double dbeta, double rtol) {
  using stan::math::var;
  var a = alpha;
  var b = beta;
  var lp = stan::math::beta_lcdf(y, a, b);
  std::vector<var> x{a, b};
  std::vector<double> g;
  lp.grad(x, g);
  EXPECT_TRUE(std::isfinite(g[0]))
      << "y=" << y << " a=" << alpha << " b=" << beta;
  EXPECT_TRUE(std::isfinite(g[1]))
      << "y=" << y << " a=" << alpha << " b=" << beta;
  EXPECT_NEAR(dalpha, g[0], rtol * std::fabs(dalpha))
      << "y=" << y << " a=" << alpha << " b=" << beta;
  EXPECT_NEAR(dbeta, g[1], rtol * std::fabs(dbeta))
      << "y=" << y << " a=" << alpha << " b=" << beta;
}

}  // namespace

TEST_F(AgradRev, ProbDistributionsBeta_lcdf_shape_gradients) {
  expect_beta_lcdf_grads(0.25, 2, 3, -0.87136898798556226, 0.40853599542368053,
                         1e-12);
  // large shapes: beta(alpha, beta) underflows
  expect_beta_lcdf_grads(0.5, 600, 600, -0.0230409429460951, 0.0230409429460951,
                         1e-12);
  // y far above the mean: the cdf is 1 and the gradients are 1e-138
  expect_beta_lcdf_grads(0.758, 352, 590, -3.3772717797123241e-138,
                         4.5390666074606647e-138, 1e-12);
  // integer beta
  expect_beta_lcdf_grads(0.5, 10, 1, -0.69314718055994531, 2.3201212683796071,
                         1e-12);
}
