#ifndef STAN_MATH_REV_FUN_IDENTITY_CONSTRAIN_HPP
#define STAN_MATH_REV_FUN_IDENTITY_CONSTRAIN_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/identity_free.hpp>
#include <stan/math/rev/fun/to_var_value.hpp>

namespace stan {
namespace math {

/**
 * Returns the result of applying the identity constraint
 * transform to the input. For any of the inputs being a `var_value` type
 * this promotes the input `x` to `var_value`.
 *
 * @tparam VarMat Any type.
 * @tparam Types Any type with one of `VarMat` and `Types` being a `var_value`
 * matrix.
 * @param[in] x object
 * @return transformed input
 */
template <typename VarMat, typename... Types,
          require_any_var_matrix_t<VarMat, Types...>* = nullptr>
inline auto identity_free(VarMat&& x, Types&&... /* args */) {
  return to_var_value(x);
}

}  // namespace math
}  // namespace stan

#endif
