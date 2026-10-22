#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <vector>

namespace {

// Shape gradients of inc_beta against a high-precision reference (mpmath,
// 80 digits), at points where a finite-difference check is not informative.
void expect_inc_beta_grads(double a, double b, double z, double dda, double ddb,
                           double rtol) {
  using stan::math::var;
  var av = a;
  var bv = b;
  var zv = z;
  var f = stan::math::inc_beta(av, bv, zv);
  std::vector<var> x{av, bv};
  std::vector<double> g;
  f.grad(x, g);
  EXPECT_TRUE(std::isfinite(g[0])) << "a=" << a << " b=" << b << " z=" << z;
  EXPECT_TRUE(std::isfinite(g[1])) << "a=" << a << " b=" << b << " z=" << z;
  EXPECT_NEAR(dda, g[0], rtol * std::fabs(dda))
      << "a=" << a << " b=" << b << " z=" << z;
  EXPECT_NEAR(ddb, g[1], rtol * std::fabs(ddb))
      << "a=" << a << " b=" << b << " z=" << z;
}

}  // namespace

TEST_F(AgradRev, inc_beta_shape_gradients) {
  expect_inc_beta_grads(2, 3, 0.25, -0.22805360232434637, 0.10692153005229139,
                        1e-12);
  // a + b above the beta(a, b) underflow limit
  expect_inc_beta_grads(600, 600, 0.5, -0.01152047147304755,
                        0.01152047147304755, 1e-12);
  // z^a (1 - z)^b underflows while I is 1
  expect_inc_beta_grads(352, 590, 0.758, -3.3772717797123242e-138,
                        4.5390666074606649e-138, 1e-12);
  // integer b
  expect_inc_beta_grads(10, 1, 0.5, -0.00067690154351557159,
                        0.00226574342615196, 1e-12);
}

TEST_F(AgradRev, inc_beta_upper_bound_infinite_density) {
  // At z = 1, I is 1 for all shapes, so the shape gradients are 0; with
  // b < 1 the density, the z gradient, is infinite.
  using stan::math::var;
  var a = 1.05;
  var b = 0.35;
  var z = 1.0;
  var f = stan::math::inc_beta(a, b, z);
  EXPECT_NO_THROW(f.grad());
  EXPECT_EQ(1.0, f.val());
  EXPECT_EQ(0.0, a.adj());
  EXPECT_EQ(0.0, b.adj());
  EXPECT_TRUE(std::isinf(z.adj()) && z.adj() > 0);
}
