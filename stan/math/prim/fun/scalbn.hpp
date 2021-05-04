#ifndef STAN_MATH_PRIM_FUN_SCALBN_HPP
#define STAN_MATH_PRIM_FUN_SCALBN_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T, typename = require_autodiff_t<T>>
double scalbn(const T& x, int n) {
  return std::scalbn(value_of_rec(x), n);
}

}  // namespace math
}  // namespace stan

#endif
