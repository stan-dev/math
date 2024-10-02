#ifndef STAN_MATH_PRIM_CONSTRAINT_IDENTITY_CONSTRAIN_HPP
#define STAN_MATH_PRIM_CONSTRAINT_IDENTITY_CONSTRAIN_HPP

#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {

/**
 * Returns the result of applying the identity constraint
 * transform to the input. This promotes the input to the least upper
 * bound of the input types.
 *
 * @tparam Jacobian if `true`, increment log density accumulator with log
 * absolute Jacobian determinant of constraining transform
 * @tparam T type of value to promote
 * @tparam Types Other types.
 * @param[in] x object
 * @return transformed input
 */
template <bool Jacobian = false, typename T, typename... Types,
          require_all_not_var_matrix_t<T, Types...>* = nullptr>
inline auto identity_constrain(T&& x, Types&&... /* args */) {
  return promote_scalar_t<return_type_t<T, Types...>, T>(x);
}

}  // namespace math
}  // namespace stan

#endif
