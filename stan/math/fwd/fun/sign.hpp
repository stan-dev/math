#ifndef STAN_MATH_FWD_FUN_SIGN_HPP
#define STAN_MATH_FWD_FUN_SIGN_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>

namespace stan {
namespace math {

template <typename T>
inline auto sign(const fvar<T>& x) {
  double z = value_of_rec(x);
  return (z == 0) ? 0 : z < 0 ? -1 : 1;
}

}  // namespace math
}  // namespace stan
#endif
