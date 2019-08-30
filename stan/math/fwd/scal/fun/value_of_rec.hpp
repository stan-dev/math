#ifndef STAN_MATH_FWD_SCAL_FUN_VALUE_OF_REC_HPP
#define STAN_MATH_FWD_SCAL_FUN_VALUE_OF_REC_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/fun/value_of_rec.hpp>
#include <type_traits>
#include <utility>

namespace stan {
namespace math {

/**
 * Return the value of the specified variable.
 *
 * T must implement value_of_rec.
 *
 * @tparam T Scalar type
 * @param x Variable.
 * @return Value of variable.
 */
template <typename T, enable_if_fvar<std::decay_t<T>>* = nullptr>
inline auto&& value_of_rec(T&& x) {
  return value_of_rec(std::forward<T>(x).val_);
}

}  // namespace math
}  // namespace stan
#endif
