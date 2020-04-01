#ifndef STAN_MATH_PRIM_FUN_GAUS_INTERP
#define STAN_MATH_PRIM_FUN_GAUS_INTERP

#include <stan/math/prim/fun/conv_gaus_line.hpp>
#include <cmath>
#include <vector>
#include <algorithm>

using namespace std;

namespace stan {
namespace math {

// set tolerance for how close interpolated values need to be to
// the specified function values
const double INTERP_TOL = 1e-8;
const double SIG2_SCALE = 0.1;

/*
given a set of points (x_i, y_i) with x_i's in increasing order, find
a_i, b_i such that the line between (x_i, y_i) and (x_{i+1}, y_{i+1}) is
f(t) = a_i*t + b_i, store the a_i's and b_i's in as and bs.
*/
void lin_interp_coefs(int n, std::vector<double> xs, std::vector<double> ys,
                      std::vector<double>& as, std::vector<double>& bs) {
  // initialize size of vectors of linear coefficients
  as.resize(n - 1);
  bs.resize(n - 1);

  // find slope and intercept between each point
  for (int i = 0; i < n - 1; i++) {
    as[i] = (ys[i + 1] - ys[i]) / (xs[i + 1] - xs[i]);
    bs[i] = -xs[i] * as[i] + ys[i];
  }
}

/*
linear interpolation at one point, x, given the coefficients
*/
double lin_interp_pt(int n, vector<double> xs, vector<double> ys,
                     vector<double> as, vector<double> bs, double x) {
  // find interval where x lives
  if (x <= xs[0])
    return ys[0];
  if (x >= xs[n - 1])
    return ys[n - 1];

  auto lb = upper_bound(xs.begin(), xs.end(), x);
  int ind = distance(xs.begin(), lb) - 1;
  ind = std::max(0, ind);

  return as[ind] * x + bs[ind];
}

/*
linear interpolation at a vector of points
*/
vector<double> lin_interp(int n, std::vector<double> xs, std::vector<double> ys,
                          int n_new, std::vector<double> xs_new) {
  // compute coefficients of linear interpolation
  vector<double> as, bs;
  lin_interp_coefs(n, xs, ys, as, bs);

  // evaluate at new points
  vector<double> ys_new;
  for (int i = 0; i < n_new; i++) {
    ys_new.push_back(lin_interp_pt(n, xs, ys, as, bs, xs_new[i]));
  }
  return ys_new;
}

/*
given a set of points and a width for the gaussian kernel, do a convolution
and evaluate at one point, x
*/
template <typename Tx>
inline return_type_t<Tx> interp_gaus_pt(int n, std::vector<double> xs,
                                        std::vector<double> ys,
                                        std::vector<double> as,
                                        std::vector<double> bs, Tx const& x,
                                        double sig2) {
  // extend out first and last lines for convolution
  double sig = std::sqrt(sig2);
  xs[0] += -10 * sig;
  xs[n - 1] += 10 * sig;

  // no need to convolve far from center of gaussian, so
  // get lower and upper indexes for integration bounds
  auto lb = upper_bound(xs.begin(), xs.end(), x - 10 * sig);
  int ind_start = distance(xs.begin(), lb) - 1;
  ind_start = std::max(0, ind_start);

  auto ub = upper_bound(xs.begin(), xs.end(), x + 10 * sig);
  int ind_end = distance(xs.begin(), ub);
  ind_end = std::min(n - 1, ind_end);

  // sum convolutions over intervals
  using return_t = return_type_t<Tx>;
  return_t y = 0;
  for (int i = ind_start; i < ind_end; i++) {
    y += conv_gaus_line(xs[i], xs[i + 1], as[i], bs[i], x, sig2);
  }

  return y;
}

/*
find the smallest difference between successive elements in a sorted vector
*/
template <typename Tx>
double min_diff(int n, std::vector<Tx> xs) {
  double dmin = value_of(xs[1]) - value_of(xs[0]);
  for (int i = 1; i < n - 1; i++) {
    if (value_of(xs[i + 1]) - value_of(xs[i]) < dmin) {
      dmin = value_of(xs[i + 1]) - value_of(xs[i]);
    }
  }
  return dmin;
}

/*
given a set of pairs (x_i, y_i), do a gaussian interpolation through those
points and evaluate the interpolation at the points xs_new
*/
template <typename Tx>
inline vector<Tx> gaus_interp(int n, std::vector<double> xs,
                              std::vector<double> ys, int n_new,
                              std::vector<Tx> xs_new) {
  // find minimum distance between points for std of gaussian kernel
  double sig2 = square(min_diff(n, xs) * SIG2_SCALE);

  // copy ys into a new vector
  std::vector<double> y2s;
  y2s.resize(n);
  for (int i = 0; i < n; i++) {
    y2s[i] = ys[i];
  }

  // interatively find interpolation that coincides with ys at xs
  int max_iters = 50;
  double dmax, dd;
  std::vector<double> as, bs;
  for (int j = 0; j < max_iters; j++) {
    // linear interpolation for new ys
    lin_interp_coefs(n, xs, y2s, as, bs);
    dmax = 0;
    for (int i = 0; i < n; i++) {
      dd = ys[i] - interp_gaus_pt(n, xs, y2s, as, bs, xs[i], sig2);
      y2s[i] += dd;
      dmax = std::max(std::abs(dd), dmax);
    }
    if (dmax < INTERP_TOL)
      break;
  }

  // fill ys_new with interpolated values at new_xs
  std::vector<Tx> ys_new;
  for (int i = 0; i < n_new; i++) {
    ys_new.push_back(interp_gaus_pt(n, xs, y2s, as, bs, xs_new[i], sig2));
  }
  return ys_new;
}

}  // namespace math
}  // namespace stan

#endif
