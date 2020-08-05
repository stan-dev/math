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

/*
find the smallest difference between successive elements in a sorted vector
*/
template <typename Tx>
double min_diff(int n, std::vector<Tx> const& xs) {
  double dmin = value_of(xs[1]) - value_of(xs[0]);
  for (int i = 1; i < n - 1; i++) {
    if (value_of(xs[i + 1]) - value_of(xs[i]) < dmin) {
      dmin = value_of(xs[i + 1]) - value_of(xs[i]);
    }
  }
  return dmin;
}

} // namespace internal

/**
 * Given a set of reference points \f$(xs_i, ys_i)\f$, create a mollifier
 * that intersects the reference points. This function requires as input
 * a struct created by the function gaus_interp_precomp. The algorithm
 * used to create the mollifier is an iterative algorithm that works 
 * as follows. First a linear 
 * interpolation is created through the reference points. Then, the 
 * linear interpolation is convolved with a Gaussian whose width is 
 * proportional the smallest distance between successive points 
 * \f$xs_i\f$ and \f$xs_{i+1}\f$. Since the convolution is unlikely to 
 * intersect the reference points, the y-values of the reference points
 * are shifted and the process repeats itself until the interpolation 
 * intersects all reference points. 
 *
 * @param xs vector of independent variable of reference points
 * @param ys vector of dependent variable of reference points
 * @param params a struct created by gaus_interp_precomp that 
 * @param x the point at which to evaluate the interpolation
 * @return value of the interpolation at x
 */
template <typename Tx>
inline return_type_t<Tx> gaus_interp(std::vector<double> const& xs,
				     std::vector<double> const& ys,
				     gaus_interp_params const& params,
				     Tx const& x) {
  const double NSTDS = 10;

  // enforce that interpolation point is between smallest and largest 
  // reference point
  static char const* function = "interp_gaus_pt";
  check_less_or_equal(function, "Interpolation point", x, xs.back());
  check_greater_or_equal(function, "Interpolation point", x, xs.front());

  int n = xs.size();

  // create copy of xs so that endpoints can be extended
  std::vector<double> xs2 = xs;

  // extend out first and last lines for convolution
  double sig = std::sqrt(params.sig2);
  xs2[0] = xs[0] - NSTDS * sig;
  xs2[n - 1] = xs[n-1] + NSTDS * sig;

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
    y += conv_gaus_line(xs2[i], xs2[i+1], params.as[i], params.bs[i], x, 
			params.sig2);
  }

  return y;
}

/**
 * This function was written to be used with gaus_interp. This function
 * computes the shifted y-values of the reference points of an interpolation
 * in such a way that when that piecewise linear function is convolved 
 * with a Gaussian kernel, the resulting function coincides with the 
 * points \f$(xs_i, ys_i)\f$ inputted into this function. The output of this
 * function depends heavily on the choice of width of the Gaussian 
 * kernel, which at the time of writing, is set to one tenth the 
 * minimum distance between successive elements of the vector xs. 
 * A tolerance for the maximum distance between the interpolation and 
 * all reference points is also set manually and is not an input. 
 *
 *
 * @param xs vector of independent variable of reference points
 * @param ys vector of dependent variable of reference points
 * @param params a struct created by gaus_interp_precomp that 
 * @param x the point at which to evaluate the interpolation
 * @return struct containing slopes, intercepts, and width of kernel 
 */
gaus_interp_params gaus_interp_precomp(std::vector<double> const& xs,
				       std::vector<double> const& ys) {
  using internal::min_diff;
  gaus_interp_params params;
  internal::asbs coefs;
  const double INTERP_TOL = 1e-8;
  const double SIG2_SCALE = 0.1;
  int n = xs.size();

  // find minimum distance between points for std of gaussian kernel
  params.sig2 = square(min_diff(n, xs) * SIG2_SCALE);
  params.as.resize(n-1);
  params.bs.resize(n-1);

  // copy ys into a new vector that will be changed
  std::vector<double> y2s = ys;

  // interatively find interpolation that coincides with ys at xs
  int max_iters = 50;
  double dmax, dd;
  for (int j = 0; j < max_iters; j++) {
    // find slope and intercept of line between each point
    for (int i = 0; i < n - 1; i++) {
      params.as[i] = (y2s[i + 1] - y2s[i]) / (xs[i + 1] - xs[i]);
      params.bs[i] = -xs[i] * params.as[i] + y2s[i];
    }

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

/**
 * This function combines gaus_interp_precomp and gaus_interp.
 * It takes as input two vectors of reference points (xs and ys)
 * in addition to a vector, xs_new, of points at which the 
 * function will evaluate the interpolation through those reference 
 * points.
 *
 * @param xs vector of independent variable of reference points
 * @param ys vector of dependent variable of reference points
 * @param xs_new vector of point at which to evaluate interpolation
 * @return vector of interpolation values
 */

template <typename Tx>
inline vector<Tx> gaus_interp_vect(std::vector<double> const& xs,
				   std::vector<double> const& ys,
				   std::vector<Tx> const& xs_new) {
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
