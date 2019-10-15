#ifndef STAN_MATH_FWD_SCAL_FUN_TO_FVAR_HPP
#define STAN_MATH_FWD_SCAL_FUN_TO_FVAR_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>

namespace stan {
namespace math {

template <typename T>
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
template <typename T>
inline const fvar<T>& to_fvar(const fvar<T>& x) {
  return x;
}

/**
 * Specialization of to_fvar for non-const fvars
 *
 *
 * @param[in,out] x A forward automatic differentation variables.
 * @return The input forward automatic differentiation variables.
 */
template <typename T>
inline fvar<T>& to_fvar(fvar<T>& x) {
  return x;
}

}  // namespace math
}  // namespace stan
#endif
