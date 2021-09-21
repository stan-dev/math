#include <test/unit/math/test_ad.hpp>
#include <limits>
#include <vector>

double ABS_TOL = 1e-12;

TEST(mathPrimInterpLin, throwing) {
  using stan::math::interp_lin;
  double nan = std::numeric_limits<double>::quiet_NaN();
  double x = 0.5;
  std::vector<double> xs, ys;

  // check that when xs are not increasing, an error throws
  int n = 2;
  xs = {1, 1};
  ys = {0, 2};
  EXPECT_THROW(interp_lin(xs, ys, x), std::domain_error);

  // check when xs contain a nan
  xs = {nan, 1};
  ys = {0, 2};
  EXPECT_THROW(interp_lin(xs, ys, x), std::domain_error);

  // xs must contain at least two elements
  xs = {1};
  ys = {0, 2};
  EXPECT_THROW(interp_lin(xs, ys, x), std::domain_error);

  // ys can't contain nan
  xs = {0, 1};
  ys = {0, nan};
  x = 0.5;
  EXPECT_THROW(interp_lin(xs, ys, x), std::domain_error);

  // x can't be nan
  xs = {0, 1};
  ys = {0, 2};
  x = nan;
  EXPECT_THROW(interp_lin(xs, ys, x), std::domain_error);
}

TEST(mathPrimInterpLin, interp_line) {
  using stan::math::interp_lin;

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

  // interpolate
  std::vector<double> ys_new(n_interp);
  for (int i = 0; i < n_interp; i++) {
    ys_new[i] = interp_lin(xs, ys, xs_new[i]);
  }

  // test points
  double tmp, y;
  for (int i = 0; i < n_interp; i++) {
    tmp = (ys[1] - ys[0]) / (xs[1] - xs[0]);
    y = tmp * xs_new[i] + ys[0] - tmp * xs[0];
    ASSERT_NEAR(ys_new[i], y, ABS_TOL);
  }

  // check values outside of range of reference points
  ASSERT_NEAR(interp_lin(xs, ys, -1), ys[0], ABS_TOL);
  ASSERT_NEAR(interp_lin(xs, ys, 100), ys[1], ABS_TOL);

  // xs with more than 2 elements
  std::vector<double> xs2 = {0, 1, 2, 3, 4, 5, 6};
  std::vector<double> ys2 = {0, 2, 5, 2, 3, 2, 2};
  double x = 0.5;
  ASSERT_NEAR(interp_lin(xs, ys, x), interp_lin(xs2, ys2, x), ABS_TOL);
}

TEST(mathPrimInterpLin, matching_reference_interp_pts) {
  using stan::math::interp_lin;

  // check that interpolation returns the same function
  // when interpolation points are the same as reference points

  // generate function tabulation
  int n = 3;
  std::vector<double> xs = {0, 1, 2};
  std::vector<double> ys = {0, 2, 1};

  // create interpolation points, same as reference points
  std::vector<double> xs_new = xs;

  // create interpolation
  std::vector<double> ys_new = interp_lin(xs, ys, xs_new);

  // test points
  for (int i = 0; i < n; i++) {
    ASSERT_NEAR(ys_new[i], ys[i], ABS_TOL);
  }
}
