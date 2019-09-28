#ifndef STAN_MATH_PRIM_SCAL_FUN_INV_SQRT_HPP
#define STAN_MATH_PRIM_SCAL_FUN_INV_SQRT_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/fun/inv.hpp>
#include <cmath>

namespace stan {
namespace math {

inline double inv_sqrt(double x) {
  using std::sqrt;
  return inv(sqrt(x));
}

}  // namespace math
}  // namespace stan

#endif
