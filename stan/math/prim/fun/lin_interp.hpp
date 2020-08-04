#ifndef STAN_MATH_PRIM_FUN_LIN_INTERP
#define STAN_MATH_PRIM_FUN_LIN_INTERP

#include <cmath>
#include <vector>
#include <algorithm>

using namespace std;

namespace stan {
namespace math {

double lin_interp(std::vector<double> xs, std::vector<double> ys, double x) {
  double x1, x2, y1, y2;
  int n;

  n = xs.size();

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
  if (x == *(ub-1)) return ys[ind-1];

  // do linear interpolation 
  x1 = *(ub-1);
  x2 = *ub;
  y1 = ys[ind-1];
  y2 = ys[ind];
  
  return y1 + (x - x1) * (y2 - y1) / (x2-x1);
}


}  // namespace math
}  // namespace stan

#endif
