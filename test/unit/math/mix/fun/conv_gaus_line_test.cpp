#include <test/unit/math/test_ad.hpp>
#include <limits>
#include <vector>

/*
check derivative using finite differences and autodiffing
*/
TEST(mathMixGausInterp, gaus_conv_der) {
  double t0, t1, a, b, x0, sig2, y, h, yn, yp, x0p, x0n, dder, dder2;
  double alpha, sig;
  using stan::math::conv_gaus_line;

  t0 = 0;
  t1 = 5;
  a = -1;
  b = 2;
  x0 = 0.99;
  sig2 = 2;

  // derivative with finite difference
  h = 1e-4;
  x0n = x0 + h;
  x0p = x0 - h;
  yn = conv_gaus_line(t0, t1, a, b, x0n, sig2);
  yp = conv_gaus_line(t0, t1, a, b, x0p, sig2);
  dder = (yn - yp) / (2 * h);

  // derivative with autodiff
  using stan::math::var;
  var v = x0;
  var vy = conv_gaus_line(t0, t1, a, b, v, sig2);
  vy.grad();

  ASSERT_NEAR(dder, v.adj(), 1e-5);
}
