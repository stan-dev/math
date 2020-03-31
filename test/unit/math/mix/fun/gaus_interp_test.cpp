#include <test/unit/math/test_ad.hpp>
#include <limits>
#include <vector>

double ABS_TOL = 1e-12;

TEST(mathMixGausInterp, interp_line) {
  using stan::math::gaus_interp;

  // check that interpolation of line returns the same function
  // generate function tabulation
  int n = 2;
  double xmin = 0;
  double xmax = 1;
  vector<double> xs = {0, 1};
  vector<double> ys = {0, 2};

  // vector of pts at which to compute interpolation
  vector<double> xs_new;
  int n_interp = 100;
  double t0 = xmin;
  double t1 = xmax;
  double t;
  for (int i=0; i < n_interp; i++) {
    t = t0 + i * (t1 - t0) / (n_interp - 1);
    xs_new.push_back(t);
  }

  // create interpolation
  vector<double> ys_new;
  ys_new = gaus_interp(n, xs, ys, n_interp, xs_new);

  // test points
  double tmp, y;
  for (int i=0; i < n_interp; i++) {
    tmp = (ys[1] - ys[0]) / (xs[1] - xs[0]);
    y = tmp * xs_new[i] + ys[0] - tmp * xs[0];
    ASSERT_NEAR(ys_new[i], y, ABS_TOL);
  }
}

TEST(mathMixGausInterp, gaus_and_lin_interp) {
  using stan::math::gaus_interp;
  using stan::math::lin_interp;

  // check that interpolation of line returns the same function
  // generate function tabulation
  int n = 30;
  double xmin = 0;
  double xmax = 1;
  double x;
  vector<double> xs, ys;
  for (int i=0; i < n; i++) {
    x = xmin + i * (xmax - xmin) / (n - 1);
    xs.push_back(x);
    ys.push_back(x*x);
  }

  // vector of pts at which to compute interpolation
  vector<double> xs_new;
  int n_interp = 100;
  double t0 = xmin;
  double t1 = xmax;
  double t;
  for (int i=0; i < n_interp; i++) {
    t = t0 + i * (t1 - t0) / (n_interp - 1);
    xs_new.push_back(t);
  }

  // create interpolation
  vector<double> ys_new_gaus, ys_new_lin;
  ys_new_gaus = gaus_interp(n, xs, ys, n_interp, xs_new);
  ys_new_lin = lin_interp(n, xs, ys, n_interp, xs_new);

  // test points
  double tmp, y;
  for (int i=0; i < n_interp; i++) {
    ASSERT_NEAR(ys_new_lin[i], ys_new_gaus[i], 1e-4);
  }
}


TEST(mathMixGausInterp, derivs) {
  using stan::math::var;
  using stan::math::gaus_interp;
  std::vector<double> xs, ys, x2s, y2s, ts, as, bs, y3s, t3s;
  double xmin, xmax, x, y, x2, y2, t, t0, t1, dd, dder, dder2;
  double x0, x1, y0, y1, tmp;
  int n;

  // generate function tabulation
  n = 10;
  xmin = 0;
  xmax = 1;
  for (int i = 0; i < n; i++) {
    x = xmin + i * (xmax - xmin) / (n - 1);
    xs.push_back(x);
    y = rand() % 100;
    ys.push_back(y);
  }

  // create vector of interpolation pts
  std::vector<stan::math::var> xs_new_v;
  std::vector<double> xs_new;
  int n_interp = 100;
  t0 = xs[1] - 1e-4;
  t0 = xs[0];
  t1 = xmax;
  for (int i=0; i < n_interp; i++) {
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
  for (int i=0; i < n_interp; i++) {
    ys_new_v[i].grad();
    dder = xs_new_v[i].adj();
    ys_new_dder.push_back(dder);
  }

  // take derivative of interpolation using finite differencing
  double h = 1e-4;
  xs_new_p = xs_new;
  xs_new_n = xs_new;
  for (int i=0; i < n_interp; i++) {
    xs_new_p[i] += -h;
    xs_new_n[i] += h;
    ys_new_p = gaus_interp(n, xs, ys, n_interp, xs_new_p);
    ys_new_n = gaus_interp(n, xs, ys, n_interp, xs_new_n);
    dder2 = (ys_new_n[i] - ys_new_p[i]) / (2 * h);
  }


}
