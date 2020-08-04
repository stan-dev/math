#ifndef STAN_MATH_PRIM_FUN_GAUS_INTERP
#define STAN_MATH_PRIM_FUN_GAUS_INTERP

#include <cmath>
#include <vector>
#include <algorithm>
#include <stan/math/prim/fun/conv_gaus_line.hpp>

using namespace std;

namespace stan {
namespace math {

// struct for passing precomputation data
struct gaus_interp_params {
  std::vector<double> as;
  std::vector<double> bs;
  double sig2;
};

namespace internal {

// struct for passing precomputation data
struct asbs {
  std::vector<double> as;
  std::vector<double> bs;
};

// set tolerance for how close interpolated values need to be to
// the specified function values
const double INTERP_TOL = 1e-8;
const double SIG2_SCALE = 0.1;
const double NSTDS = 10;

/*
given a set of points (x_i, y_i) with x_i's in increasing order, find
a_i, b_i such that the line between (x_i, y_i) and (x_{i+1}, y_{i+1}) is
f(t) = a_i*t + b_i, store the a_i's and b_i's in as and bs.
*/
asbs lin_interp_coefs(std::vector<double> xs, std::vector<double> ys) {
  
  int n = xs.size();

  // initialize size of vectors of linear coefficients
  std::vector<double> as(n - 1);
  std::vector<double> bs(n - 1);

  // find slope and intercept between each point
  for (int i = 0; i < n - 1; i++) {
    as[i] = (ys[i + 1] - ys[i]) / (xs[i + 1] - xs[i]);
    bs[i] = -xs[i] * as[i] + ys[i];
  }
  // return struct with as and bs
  asbs coefs;
  coefs.as = as;
  coefs.bs = bs;
  return coefs;
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

} // namespace internal

/*
given a set of points and a width for the gaussian kernel, do a convolution
and evaluate at one point, x
*/
template <typename Tx>
inline return_type_t<Tx> gaus_interp(std::vector<double> xs,
				     std::vector<double> ys,
				     gaus_interp_params params,
				     Tx const& x) {
  using internal::NSTDS;

  // enforce that interpolation point is between smallest and largest 
  // reference point
  static char const* function = "interp_gaus_pt";
  check_less_or_equal(function, "Interpolation point", x, xs.back());
  check_greater_or_equal(function, "Interpolation point", x, xs.front());

  int n = xs.size();

  // extend out first and last lines for convolution
  double sig = std::sqrt(params.sig2);
  xs[0] += -NSTDS * sig;
  xs[n - 1] += NSTDS * sig;

  // no need to convolve far from center of gaussian, so
  // get lower and upper indexes for integration bounds
  auto lb = upper_bound(xs.begin(), xs.end(), x - NSTDS * sig);
  int ind_start = distance(xs.begin(), lb) - 1;
  ind_start = std::max(0, ind_start);

  auto ub = upper_bound(xs.begin(), xs.end(), x + NSTDS * sig);
  int ind_end = distance(xs.begin(), ub);
  ind_end = std::min(n - 1, ind_end);

  // sum convolutions over intervals
  using return_t = return_type_t<Tx>;
  return_t y = 0;
  for (int i = ind_start; i < ind_end; i++) {
    y += conv_gaus_line(xs[i], xs[i+1], params.as[i], params.bs[i], x, 
			params.sig2);
  }

  return y;
}

/*
given a set of pairs (x_i, y_i), do a gaussian interpolation through those
points and evaluate the interpolation at the points xs_new
*/
gaus_interp_params gaus_interp_precomp(std::vector<double> xs,
				       std::vector<double> ys) {
  using internal::min_diff;
  using internal::lin_interp_coefs;
  using internal::SIG2_SCALE;
  using internal::INTERP_TOL;
  gaus_interp_params params;
  internal::asbs coefs;
  int n = xs.size();

  // find minimum distance between points for std of gaussian kernel
  params.sig2 = square(min_diff(n, xs) * SIG2_SCALE);

  // copy ys into a new vector
  std::vector<double> y2s;
  y2s.resize(n);
  for (int i = 0; i < n; i++) {
    y2s[i] = ys[i];
  }

  // interatively find interpolation that coincides with ys at xs
  int max_iters = 50;
  double dmax, dd;
  for (int j = 0; j < max_iters; j++) {
    // linear interpolation for new ys
    coefs = lin_interp_coefs(xs, y2s);
    params.as = coefs.as;
    params.bs = coefs.bs;
    dmax = 0;
    for (int i = 0; i < n; i++) {
      dd = ys[i] - gaus_interp(xs, y2s, params, xs[i]);
      y2s[i] += dd;
      dmax = std::max(std::abs(dd), dmax);
    }
    if (dmax < INTERP_TOL)
      break;
  }

  return params;
}
template <typename Tx>
inline vector<Tx> gaus_interp_vect(std::vector<double> xs,
				   std::vector<double> ys,
				   std::vector<Tx> xs_new) {
  int n_interp = xs_new.size();
  std::vector<Tx> ys_new(n_interp);

  // create interpolation
  gaus_interp_params params = gaus_interp_precomp(xs, ys);
  for (int i=0; i<n_interp; i++) {
    ys_new[i] = gaus_interp(xs, ys, params, xs_new[i]);
  }
  return ys_new;
}

}  // namespace math
}  // namespace stan

#endif
