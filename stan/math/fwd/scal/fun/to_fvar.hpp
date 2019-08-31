#ifndef STAN_MATH_FWD_SCAL_FUN_TO_FVAR_HPP
#define STAN_MATH_FWD_SCAL_FUN_TO_FVAR_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/meta.hpp>
#include <utility>

namespace stan {
namespace math {

template <typename T,
          std::enable_if_t<is_var<std::decay_t<T>>::value
                           || std::is_arithmetic<std::decay_t<T>>::value>...>
inline fvar<T> to_fvar(const T& x) {
  return fvar<T>(x);
}

/**
 * Specialization of to_fvar for const fvars
 *
 *
 * @param[in,out] x A forward automatic differentation variables.
 * @return The input forward automatic differentiation variables.
 */
template <typename T, fvar_type<T>...>
inline const auto& to_fvar(T&& x) {
  return std::forward<T>(x);
}

}  // namespace math
}  // namespace stan
#endif
