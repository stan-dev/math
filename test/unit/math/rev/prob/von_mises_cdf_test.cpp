#include <stan/math/rev.hpp>
#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <vector>
#include <cmath>

double ABS_TOL = 1e-12;

TEST(VonMises, von_mises_gradients) {
  double y;
  stan::math::var a = 4.0;
  stan::math::var b = stan::math::sqrt(a);
  b.grad();

  using stan::math::von_mises_cdf;
  using stan::math::von_mises_lpdf;

  stan::math::var x = 0.5;
  stan::math::var mu = 0.0;
  stan::math::var k = 1.0;
  stan::math::var f = von_mises_cdf(x, mu, k);
  f.grad();

  // compute partial derivatives analytically for mu, x
  y = von_mises_lpdf(value_of(x), value_of(mu), value_of(k));
  y = exp(y);
  ASSERT_NEAR(x.adj(), y, ABS_TOL);
  ASSERT_NEAR(mu.adj(), -y, ABS_TOL);

  // compute derivative using finite differences
  double h = 1e-5;
  double xd = value_of(x);
  double mud = value_of(mu);
  double kd = value_of(k);
  double fn = von_mises_cdf(xd, mud, kd+h);
  double fp = von_mises_cdf(xd, mud, kd-h);
  double dder = (fn - fp) / (2 * h);
  ASSERT_NEAR(k.adj(), dder, 1e-8);

  // test another value of k, that uses the normal approximation
  k = 20.0;
  f = von_mises_cdf(x, mu, k);
  f.grad();
  kd = value_of(k);
  fn = von_mises_cdf(xd, mud, kd+h);
  fp = von_mises_cdf(xd, mud, kd-h);
  dder = (fn - fp) / (2 * h);
  ASSERT_NEAR(k.adj(), dder, 1e-8);
}
