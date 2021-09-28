#ifndef STAN_MATH_PRIM_FUN_IDENTITY_FREE_HPP
#define STAN_MATH_PRIM_FUN_IDENTITY_FREE_HPP

#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {

/**
 * Returns the result of applying the inverse of the identity
 * constraint transform to the input. This promotes the input to the least upper
 * bound of the input types.
 *
 *
 * @tparam T type of value to promote
 * @tparam Types Other types.
 * @param[in] x value to promote
 * @return value
 */
template <typename T, typename... Types,
          require_all_not_var_matrix_t<T, Types...>* = nullptr>
inline auto identity_free(T&& x, Types&&... /* args */) {
  return promote_scalar_t<return_type_t<T, Types...>, T>(x);
}

}  // namespace math
}  // namespace stan

#endif
