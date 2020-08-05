#include <test/unit/math/test_ad.hpp>
#include <limits>
#include <vector>

double ABS_TOL = 1e-12;

TEST(mathPrimGausInterp, interp_line) {
  using stan::math::gaus_interp;
  using stan::math::gaus_interp_params;
  using stan::math::gaus_interp_precomp;
  using stan::math::gaus_interp_vect;

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
  for (int i = 0; i < n_interp; i++) {
    t = t0 + i * (t1 - t0) / (n_interp - 1);
    xs_new.push_back(t);
  }

  // create interpolation using precomp
  vector<double> ys_new_gaus(n_interp);
  gaus_interp_params params = gaus_interp_precomp(xs, ys);
  for (int i = 0; i < n_interp; i++) {
    ys_new_gaus[i] = gaus_interp(xs, ys, params, xs_new[i]);
  }

  // create interpolation without precomp
  vector<double> ys_new_gaus2 = gaus_interp_vect(xs, ys, xs_new);

  // test points
  double tmp, y;
  for (int i = 0; i < n_interp; i++) {
    tmp = (ys[1] - ys[0]) / (xs[1] - xs[0]);
    y = tmp * xs_new[i] + ys[0] - tmp * xs[0];
    ASSERT_NEAR(ys_new_gaus[i], y, ABS_TOL);
    ASSERT_NEAR(ys_new_gaus2[i], y, ABS_TOL);
  }
}

TEST(mathPrimGausInterp, gaus_and_lin_interp) {
  using stan::math::gaus_interp;
  using stan::math::gaus_interp_params;
  using stan::math::gaus_interp_precomp;
  using stan::math::lin_interp;

  // check that interpolation of line returns the same function
  // generate function tabulation
  int n = 30;
  double xmin = 0;
  double xmax = 1;
  double x;
  vector<double> xs, ys;
  for (int i = 0; i < n; i++) {
    x = xmin + i * (xmax - xmin) / (n - 1);
    xs.push_back(x);
    ys.push_back(x * x);
  }

  // vector of pts at which to compute interpolation
  vector<double> xs_new;
  int n_interp = 100;
  double t0 = xmin;
  double t1 = xmax;
  double t;
  for (int i = 0; i < n_interp; i++) {
    t = t0 + i * (t1 - t0) / (n_interp - 1);
    xs_new.push_back(t);
  }

  // create interpolation
  vector<double> ys_new_gaus(n_interp);
  vector<double> ys_new_lin(n_interp);

  // linear interpolation
  ys_new_lin.resize(n_interp);
  for (int i = 0; i < n_interp; i++) {
    ys_new_lin[i] = lin_interp(xs, ys, xs_new[i]);
  }

  // gaus interpolation
  gaus_interp_params params = gaus_interp_precomp(xs, ys);
  for (int i = 0; i < n_interp; i++) {
    ys_new_gaus[i] = gaus_interp(xs, ys, params, xs_new[i]);
  }

  // test points
  double tmp, y;
  for (int i = 0; i < n_interp; i++) {
    ASSERT_NEAR(ys_new_lin[i], ys_new_gaus[i], 1e-4);
  }
}
