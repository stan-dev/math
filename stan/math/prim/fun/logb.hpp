#ifndef STAN_MATH_PRIM_FUN_LOGB_HPP
#define STAN_MATH_PRIM_FUN_LOGB_HPP

#include <stan/math/prim/meta.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T, typename = require_autodiff_t<T>>
double logb(const T& x) {
  return std::logb(value_of_rec(x));
}

}  // namespace math
}  // namespace stan

#endif
