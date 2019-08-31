#ifndef STAN_MATH_REV_SCAL_FUN_TO_VAR_HPP
#define STAN_MATH_REV_SCAL_FUN_TO_VAR_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <type_traits>
#include <utility>

namespace stan {
namespace math {

/**
 * Converts argument to an automatic differentiation variable.
 *
 * Returns a var variable with the input value.
 *
 * @param[in] x A scalar value
 * @return An automatic differentiation variable with the input value.
 */
inline var to_var(double x) { return var(x); }

/**
 * Specialization of to_var for non-const var input
 *
 *
 * @param[in,out] x An automatic differentiation variable.
 * @return The input automatic differentiation variable.
 */
template <typename T, require_var<std::decay_t<T>>...>
inline auto&& to_var(T&& x) {
  return std::forward<T>(x);
}

}  // namespace math
}  // namespace stan
#endif
