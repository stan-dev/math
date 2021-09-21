#include <test/unit/math/test_ad.hpp>

TEST(mathMixGausInterp, derivs) {
  using stan::math::interp_gauss;
  using stan::math::var;
  std::vector<double> xs, ys, ts;
  double xmin, xmax, x, y, t, t0, t1, dder, dder2, x0, x1;
  int n;

  // generate function tabulation
  n = 5;
  xmin = 0;
  xmax = 1;
  xs.resize(n);
  ys.resize(n);
  for (int i = 0; i < n; i++) {
    x = xmin + i * (xmax - xmin) / (n - 1);
    xs[i] = x;
  }
  ys[0] = 0;
  ys[1] = 0;
  ys[2] = 1;
  ys[3] = 0;
  ys[4] = 0.01;

  // create vector of interpolation pts
  std::vector<stan::math::var> xs_new_v;
  std::vector<double> xs_new;
  int n_interp = 100;

  // add a cushion to each endpoint so that we don't leave interval
  // when taking finite differences
  t0 = xs[0] + 0.01;
  t1 = xmax - 0.01;
  for (int i = 0; i < n_interp; i++) {
    t = t0 + i * (t1 - t0) / (n_interp - 1);
    xs_new.push_back(t);
    stan::math::var tvar = t;
    xs_new_v.push_back(tvar);
  }

  std::vector<double> ys_new, ys_new_p, ys_new_n, xs_new_p, xs_new_n;
  ys_new = interp_gauss(xs, ys, xs_new);

  // autodiff at each interpolation pt
  std::vector<stan::math::var> ys_new_v;
  ys_new_v = interp_gauss(xs, ys, xs_new_v);

  std::vector<double> ys_new_dder;
  for (int i = 0; i < n_interp; i++) {
    ys_new_v[i].grad();
    dder = xs_new_v[i].adj();
    ys_new_dder.push_back(dder);
  }

  // take derivative of interpolation using finite differencing
  double h = 1e-6;
  xs_new_p = xs_new;
  xs_new_n = xs_new;
  for (int i = 0; i < n_interp; i++) {
    xs_new_p[i] += -h;
    xs_new_n[i] += h;
    ys_new_p = interp_gauss(xs, ys, xs_new_p);
    ys_new_n = interp_gauss(xs, ys, xs_new_n);
    dder2 = (ys_new_n[i] - ys_new_p[i]) / (2 * h);
    ASSERT_NEAR(ys_new_dder[i], dder2, 1e-5);
  }
}
