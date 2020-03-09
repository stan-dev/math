#ifndef STAN_MATH_PRIM_FUN_IDENTITY_FREE_HPP
#define STAN_MATH_PRIM_FUN_IDENTITY_FREE_HPP

#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {

/**
 * Returns the result of applying the inverse of the identity
 * constraint transform to the input.
 *
 * <p>This function is a no-op and mainly useful as a placeholder
 * in auto-generated code.
 *
 * @tparam Scalar type of value
 * @param[in] x value
 * @return value
 */
template <typename Scalar, require_all_stan_scalar_t<Scalar>* = nullptr>
inline decltype(auto) identity_free(Scalar&& x) {
  return std::forward<Scalar>(x);
}

/**
 * Returns the result of applying the inverse of the identity
 * constraint transform to the input.
 *
 * This function is for N-ary op functions to make sure type promotion happens
 * correctly, thus allowing the use of auto.
 *
 *
 * @tparam Scalar type of value
 * @param[in] x value
 * @return value
 */
template <typename Scalar, typename... Types,
 require_all_stan_scalar_t<Scalar, Types...>* = nullptr>
inline auto identity_free(Scalar&& x, Types&&... args) {
  return return_type_t<Scalar, Types...>(x);
}

}  // namespace math
}  // namespace stan

#endif
