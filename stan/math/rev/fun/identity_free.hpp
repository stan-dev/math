#ifndef STAN_MATH_REV_FUN_IDENTITY_CONSTRAIN_HPP
#define STAN_MATH_REV_FUN_IDENTITY_CONSTRAIN_HPP

#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {

/**
 * Returns the result of applying the identity constraint
 * transform to the input.
 *
 * <p>This method is effectively a no-op and is mainly useful as a
 * placeholder in auto-generated code.
 *
 * @tparam T Any type.
 * @param[in] x object
 * @return transformed input
 */
template <typename T, typename... Types, require_any_var_matrix_t<T, Types...>* = nullptr>
inline auto identity_free(T&& x, Types&&... args) {
  return to_var_value_if<is_any_var_matrix<T, Types...>::value>(x);
}

}  // namespace math
}  // namespace stan

#endif
