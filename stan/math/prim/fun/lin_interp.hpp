#ifndef STAN_MATH_PRIM_FUN_LIN_INTERP
#define STAN_MATH_PRIM_FUN_LIN_INTERP

#include <stan/math/prim/err.hpp>
#include <cmath>
#include <vector>
#include <algorithm>

namespace stan {
namespace math {

/**
 * This function performs linear interpolation. The function takes as
 * input two vectors of reference points \f$(xs_i, ys_i)\f$ in addition
 * to one value, \f$x\f$, at which the value of the linear interpolation
 * is to be returned. The vector xs is required to be in increasing order
 * For values of \f$x\f$ less than xs[0], the function returns ys[0].
 * For values of \f$x\f$ greater than xs[n-1], the function returns
 * ys[n-1], where xs is of length n.
 *
 * @param xs vector of independent variable of reference points
 * @param ys vector of dependent variable of reference points
 * @param x the point at which to evaluate the interpolation
 * @return value of linear interpolation at x
 */
inline double lin_interp(const std::vector<double>& xs,
                         const std::vector<double>& ys, double x) {
  static char const* function = "lin_interp";
  check_ordered(function, "xs", xs);
  check_not_nan(function, "xs", xs);
  check_not_nan(function, "ys", ys);
  check_not_nan(function, "x", x);
  check_greater(function, "xs", xs.size(), 1);
  int n = xs.size();

  // if x is less than left endpoint or greater than right, return endpoint
  if (x <= xs[0]) {
    return ys[0];
  }
  if (x >= xs[n - 1]) {
    return ys[n - 1];
  }

  // find in between which points the input, x, lives
  auto ub = std::upper_bound(xs.begin(), xs.end(), x);
  auto ind = std::distance(xs.begin(), ub);

  // check if the interpolation point falls on a reference point
  if (x == xs[ind - 1]) {
    return ys[ind - 1];
  }

  // do linear interpolation
  double x1 = xs[ind-1];
  double x2 = xs[ind];
  double y1 = ys[ind - 1];
  double y2 = ys[ind];

  return y1 + (x - x1) * (y2 - y1) / (x2 - x1);
}

/**
 * This function takes as input two vectors of reference points (xs and ys)
 * in addition to a vector, xs_new, of points at which this
 * function will evaluate a linear interpolation through those reference
 * points.
 *
 * @param xs vector of independent variable of reference points
 * @param ys vector of dependent variable of reference points
 * @param xs_new vector of point at which to evaluate interpolation
 * @return vector of interpolation values
 */

template <typename Tx>
inline std::vector<Tx> lin_interp(const std::vector<double>& xs,
				  const std::vector<double>& ys,
				  const std::vector<Tx>& xs_new) {
  int n_interp = xs_new.size();
  std::vector<Tx> ys_new(n_interp);

  // create interpolation
  for (int i = 0; i < n_interp; i++) {
    ys_new[i] = lin_interp(xs, ys, xs_new[i]);
  }
  return ys_new;
}

}  // namespace math
}  // namespace stan

#endif
