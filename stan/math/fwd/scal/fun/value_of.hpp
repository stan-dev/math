#ifndef STAN_MATH_FWD_SCAL_FUN_VALUE_OF_HPP
#define STAN_MATH_FWD_SCAL_FUN_VALUE_OF_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/fun/value_of.hpp>
#include <type_traits>
#include <utility>

namespace stan {
namespace math {

/**
 * Return the value of the specified variable.
 *
 * @param x Variable.
 * @return Value of variable.
 */
template <typename T, enable_if_fvar<std::decay_t<T>>* = nullptr>
inline auto&& value_of(T&& x) {
  return std::forward<T>(x).val_;
}

}  // namespace math
}  // namespace stan
#endif
