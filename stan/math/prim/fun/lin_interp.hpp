#ifndef STAN_MATH_PRIM_FUN_LIN_INTERP
#define STAN_MATH_PRIM_FUN_LIN_INTERP

#include <cmath>
#include <vector>
#include <algorithm>

using namespace std;

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
double lin_interp(std::vector<double> const& xs, std::vector<double> const& ys,
                  double x) {
  int n = xs.size();

  // if x is less than left endpoint or greater than right, return endpoint
  if (x <= xs[0]) {
    return ys[0];
  }
  if (x >= xs[n - 1]) {
    return ys[n - 1];
  }

  // find in between which points the input, x, lives
  auto ub = upper_bound(xs.begin(), xs.end(), x);
  int ind = distance(xs.begin(), ub);

  // check if the interpolation point falls on a reference point
  if (x == *(ub - 1))
    return ys[ind - 1];

  // do linear interpolation
  double x1 = *(ub - 1);
  double x2 = *ub;
  double y1 = ys[ind - 1];
  double y2 = ys[ind];

  return y1 + (x - x1) * (y2 - y1) / (x2 - x1);
}

}  // namespace math
}  // namespace stan

#endif
