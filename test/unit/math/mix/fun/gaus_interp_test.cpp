#include <test/unit/math/test_ad.hpp>
#include <limits>
#include <vector>

TEST(mathMixGausInterp, derivs) {
  using stan::math::gaus_interp;
  using stan::math::var;
  std::vector<double> xs, ys, x2s, y2s, ts, as, bs, y3s, t3s;
  double xmin, xmax, x, y, x2, y2, t, t0, t1, dd, dder, dder2;
  double x0, x1, y0, y1, tmp;
  int n;
  unsigned int seed = 1;

  // generate function tabulation
  n = 10;
  xmin = 0;
  xmax = 1;
  for (int i = 0; i < n; i++) {
    x = xmin + i * (xmax - xmin) / (n - 1);
    xs.push_back(x);
    y = rand_r(&seed) % 100;
    ys.push_back(y);
  }

  // create vector of interpolation pts
  std::vector<stan::math::var> xs_new_v;
  std::vector<double> xs_new;
  int n_interp = 100;
  t0 = xs[0];
  t1 = xmax;
  for (int i = 0; i < n_interp; i++) {
    t = t0 + i * (t1 - t0) / (n_interp - 1);
    xs_new.push_back(t);
    stan::math::var tvar = t;
    xs_new_v.push_back(tvar);
  }

  std::vector<double> ys_new, ys_new_p, ys_new_n, xs_new_p, xs_new_n;
  ys_new = gaus_interp(n, xs, ys, n_interp, xs_new);

  // autodiff at each interpolation pt
  std::vector<stan::math::var> ys_new_v;
  ys_new_v = gaus_interp(n, xs, ys, n_interp, xs_new_v);

  std::vector<double> ys_new_dder;
  std::vector<var> ys_new_v2;
  for (int i = 0; i < n_interp; i++) {
    ys_new_v[i].grad();
    dder = xs_new_v[i].adj();
    ys_new_dder.push_back(dder);
  }

  // take derivative of interpolation using finite differencing
  double h = 1e-8;
  xs_new_p = xs_new;
  xs_new_n = xs_new;
  for (int i = 0; i < n_interp; i++) {
    xs_new_p[i] += -h;
    xs_new_n[i] += h;
    ys_new_p = gaus_interp(n, xs, ys, n_interp, xs_new_p);
    ys_new_n = gaus_interp(n, xs, ys, n_interp, xs_new_n);
    dder2 = (ys_new_n[i] - ys_new_p[i]) / (2 * h);
    ASSERT_NEAR(ys_new_dder[i], dder2, 1e-5);
  }
}
