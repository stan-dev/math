#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <vector>

TEST_F(AgradRev, ProbDistributionsBeta_cdf_shape_gradients) {
  // For one observation the shape partials of beta_cdf are d/da I and
  // d/db I. Reference: mpmath, 80 digits. The point has b >> a with a
  // small y.
  using stan::math::var;
  var a = 1.145;
  var b = 6786;
  var p = stan::math::beta_cdf(1.8e-4, a, b);
  std::vector<var> x{a, b};
  std::vector<double> g;
  p.grad(x, g);
  EXPECT_NEAR(0.64944133941212122, p.val(), 1e-12);
  EXPECT_NEAR(-0.38621356822721375, g[0], 1e-10);
  EXPECT_NEAR(5.8439109664263782e-5, g[1], 1e-10 * 5.8e-5);
}

TEST_F(AgradRev, ProbDistributionsBeta_cdf_mixed_shape_types) {
  // one autodiff shape with the other constant must give the same partial
  using stan::math::var;
  {
    var a = 1.145;
    var p = stan::math::beta_cdf(1.8e-4, a, 6786.0);
    p.grad();
    EXPECT_NEAR(-0.38621356822721375, a.adj(), 1e-10);
  }
  {
    var b = 6786;
    var p = stan::math::beta_cdf(1.8e-4, 1.145, b);
    p.grad();
    EXPECT_NEAR(5.8439109664263782e-5, b.adj(), 1e-10 * 5.8e-5);
  }
}
