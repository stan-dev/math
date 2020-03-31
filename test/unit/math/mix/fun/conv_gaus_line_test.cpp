#include <test/unit/math/test_ad.hpp>
#include <limits>
#include <vector>

/*
evaluate the convolution using trapezoid rule
*/
double check_int_dumb(double t0, double t1, double a, double b, double x0,
		      double sig2) {
  double x, y, t, tot;
  vector<double> ys, ts;
  double pi = stan::math::pi();
  int n;

  // evaluate function at equispaced nodes
  n = 100000;
  for (int i = 0; i < n; i++) {
    t = t0 + i * (t1 - t0) / (n - 1);
    ts.push_back(t);

    y = a * t + b;
    y = y * exp(-pow(t - x0, 2) / (2 * sig2));
    ys.push_back(y);
  }

  // sum up function evaluations
  tot = 0;
  for (int i = 0; i < n; i++) {
    tot += ys[i];
  }
  tot *= (t1 - t0) / n;
  tot /= sqrt(2 * pi * sig2);

  return tot;
}

/*
check that integral using formula is close to trapezoid approximation
*/
TEST(mathMixConvGausLin, trap_test) {
  double t0, t1, a, b, x0, sig2, y, y2;
  using stan::math::conv_gaus_line;
  t0 = 0;
  t1 = 1;
  a = 1;
  b = 0;
  x0 = 1;
  sig2 = 3;
  y = check_int_dumb(t0, t1, a, b, x0, sig2);
  y2 = conv_gaus_line(t0, t1, a, b, x0, sig2);

  ASSERT_NEAR(y, y2, 1e-5);
}

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
  dder = (yn - yp) / (2*h);

  // derivative with autodiff
  using stan::math::var;
  var v = x0;
  var vy = conv_gaus_line(t0, t1, a, b, v, sig2);
  vy.grad();

  ASSERT_NEAR(dder, v.adj(), 1e-5);
}
