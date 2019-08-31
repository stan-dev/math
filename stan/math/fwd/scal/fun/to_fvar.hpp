#ifndef STAN_MATH_FWD_SCAL_FUN_TO_FVAR_HPP
#define STAN_MATH_FWD_SCAL_FUN_TO_FVAR_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/meta.hpp>
#include <utility>

namespace stan {
namespace math {

/**
 * Convert var or arithmetic to fvar.
 *
 * @tparam T The type to become the partial type of the fvar.
 * @param[in,out] x An automatic differentation variable or arithmetic type.
 * @return The input forward automatic differentiation variables.
 */
template <typename T, require_var_or_arithmetic<T>...>
inline fvar<T> to_fvar(const T& x) {
  return fvar<T>(x);
}

/**
 * Specialization of to_fvar for input fvars
 *
 *
 * @param[in,out] x A forward automatic differentation variable.
 * @return The input forward automatic differentiation variable.
 */
template <typename T, require_fvar<T>...>
inline auto&& to_fvar(T&& x) {
  return std::forward<T>(x);
}

}  // namespace math
}  // namespace stan
#endif
