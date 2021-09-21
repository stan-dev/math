#include <test/unit/math/test_ad.hpp>
#include <limits>
#include <vector>

double ABS_TOL = 1e-12;

TEST(mathPrimInterpGauss, throwing) {
  using stan::math::interp_gauss;
  using stan::math::interp_gauss_precomp;
  double nan = std::numeric_limits<double>::quiet_NaN();
  double x;
  std::vector<double> xs, ys;

  // check that when xs are not increasing, an error throws
  int n = 2;
  xs = {1, 1};
  ys = {0, 2};
  EXPECT_THROW(interp_gauss_precomp(xs, ys), std::domain_error);

  // check that when xs contain a nan
  xs = {nan, 1};
  ys = {0, 2};
  EXPECT_THROW(interp_gauss_precomp(xs, ys), std::domain_error);

  // xs must contain at least two elements
  xs = {1};
  ys = {0, 2};
  EXPECT_THROW(interp_gauss_precomp(xs, ys), std::domain_error);

  // check that error throws when trying to interpolate out of range or nan
  xs = {0, 1};
  ys = {0, 2};
  std::vector<double> params = interp_gauss_precomp(xs, ys);
  x = 1.1;
  EXPECT_THROW(interp_gauss(xs, ys, params, x), std::domain_error);
  x = -0.1;
  EXPECT_THROW(interp_gauss(xs, ys, params, x), std::domain_error);
  x = nan;
  EXPECT_THROW(interp_gauss(xs, ys, params, x), std::domain_error);

  // ys can't contain nan
  xs = {0, 1};
  ys = {0, nan};
  x = 0.5;
  EXPECT_THROW(interp_gauss(xs, ys, params, x), std::domain_error);

  // xs can't contain nan
  xs = {0, nan};
  ys = {0, 2};
  x = 0.5;
  EXPECT_THROW(interp_gauss(xs, ys, params, x), std::domain_error);

  // xs must be increasing
  xs = {1, 1};
  ys = {0, 2};
  EXPECT_THROW(interp_gauss(xs, ys, params, x), std::domain_error);

  // xs must contain at least two elements
  xs = {1};
  ys = {0, 2};
  EXPECT_THROW(interp_gauss(xs, ys, params, x), std::domain_error);
}

TEST(mathPrimInterpGauss, interp_line) {
  using stan::math::interp_gauss;
  using stan::math::interp_gauss_precomp;

  // check that interpolation of line returns the same function
  // generate function tabulation
  int n = 2;
  double xmin = 0;
  double xmax = 1;
  std::vector<double> xs = {0, 1};
  std::vector<double> ys = {0, 2};

  // vector of pts at which to compute interpolation
  std::vector<double> xs_new;
  int n_interp = 100;
  double t0 = xmin;
  double t1 = xmax;
  double t;
  for (int i = 0; i < n_interp; i++) {
    t = t0 + i * (t1 - t0) / (n_interp - 1);
    xs_new.push_back(t);
  }

  // create interpolation using precomp
  std::vector<double> ys_new_gaus(n_interp);
  std::vector<double> params = interp_gauss_precomp(xs, ys);
  for (int i = 0; i < n_interp; i++) {
    ys_new_gaus[i] = interp_gauss(xs, ys, params, xs_new[i]);
  }

  // create interpolation without precomp
  std::vector<double> ys_new_gaus2 = interp_gauss(xs, ys, xs_new);

  // test points
  double tmp, y;
  for (int i = 0; i < n_interp; i++) {
    tmp = (ys[1] - ys[0]) / (xs[1] - xs[0]);
    y = tmp * xs_new[i] + ys[0] - tmp * xs[0];
    ASSERT_NEAR(ys_new_gaus[i], y, ABS_TOL);
    ASSERT_NEAR(ys_new_gaus2[i], y, ABS_TOL);
  }
}

TEST(mathPrimInterpGauss, matching_reference_interp_pts) {
  using stan::math::interp_gauss;

  // check that interpolation returns the same function
  // when interpolation points are the same as reference points

  // generate function tabulation
  int n = 3;
  std::vector<double> xs = {0, 1, 2};
  std::vector<double> ys = {0, 2, 1};

  // create interpolation points
  std::vector<double> xs_new = xs;

  // create interpolation
  std::vector<double> ys_new = interp_gauss(xs, ys, xs_new);

  // test points
  for (int i = 0; i < n; i++) {
    ASSERT_NEAR(ys_new[i], ys[i], 1e-8);
  }
}

TEST(mathPrimInterpGauss, interp_gauss_and_lin) {
  using stan::math::interp_gauss;
  using stan::math::interp_gauss_precomp;
  using stan::math::interp_lin;

  // check that interpolation of line returns the same function
  // generate function tabulation
  int n = 30;
  double xmin = 0;
  double xmax = 1;
  double x;
  std::vector<double> xs, ys;
  for (int i = 0; i < n; i++) {
    x = xmin + i * (xmax - xmin) / (n - 1);
    xs.push_back(x);
    ys.push_back(x * x);
  }

  // vector of pts at which to compute interpolation
  std::vector<double> xs_new;
  int n_interp = 100;
  double t0 = xmin;
  double t1 = xmax;
  double t;
  for (int i = 0; i < n_interp; i++) {
    t = t0 + i * (t1 - t0) / (n_interp - 1);
    xs_new.push_back(t);
  }

  // create interpolation
  std::vector<double> ys_new_gaus(n_interp);
  std::vector<double> ys_new_lin(n_interp);

  // linear interpolation
  ys_new_lin.resize(n_interp);
  for (int i = 0; i < n_interp; i++) {
    ys_new_lin[i] = interp_lin(xs, ys, xs_new[i]);
  }

  // gaus interpolation
  std::vector<double> params = interp_gauss_precomp(xs, ys);
  for (int i = 0; i < n_interp; i++) {
    ys_new_gaus[i] = interp_gauss(xs, ys, params, xs_new[i]);
  }

  // test points
  double tmp, y;
  for (int i = 0; i < n_interp; i++) {
    ASSERT_NEAR(ys_new_lin[i], ys_new_gaus[i], 1e-4);
  }
}
