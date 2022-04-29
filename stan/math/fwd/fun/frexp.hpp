#ifndef STAN_MATH_FWD_FUN_FREXP_HPP
#define STAN_MATH_FWD_FUN_FREXP_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>

namespace stan {
namespace math {

template <typename T>
inline auto frexp(const fvar<T>& x, int* exponent) noexcept {
  return std::frexp(value_of_rec(x), exponent);
}
}  // namespace math
}  // namespace stan
#endif
