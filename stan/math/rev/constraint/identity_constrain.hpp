#ifndef STAN_MATH_REV_CONSTRAINT_IDENTITY_CONSTRAIN_HPP
#define STAN_MATH_REV_CONSTRAINT_IDENTITY_CONSTRAIN_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/constraint/identity_constrain.hpp>
#include <stan/math/rev/fun/to_var_value.hpp>
namespace stan {
namespace math {

/**
 * Returns the result of applying the identity constraint
 * transform to the input. This specialization handles convert Eigen matrices
 * of doubles to var matrix types.
 *
 * @tparam T Any type.
 * @tparam Types Any type with one of `T` and `Types` being a `var_value`
 * matrix.
 * @param[in] x object
 * @return transformed input
 */
template <typename T, typename... Types,
          require_eigen_vt<std::is_arithmetic, T>* = nullptr,
          require_any_var_matrix_t<T, Types...>* = nullptr>
inline auto identity_constrain(T&& x, Types&&... /* args */) {
  return var_value<plain_type_t<T>>(x);
}

template <typename T, typename... Types, require_eigen_vt<is_var, T>* = nullptr,
          require_any_var_matrix_t<T, Types...>* = nullptr>
inline auto identity_constrain(T&& x, Types&&... /* args */) {
  return to_var_value(x);
}

template <typename T, typename... Types, require_var_matrix_t<T>* = nullptr,
          require_any_var_matrix_t<T, Types...>* = nullptr>
inline auto identity_constrain(T&& x, Types&&... /* args */) {
  return x;
}

}  // namespace math
}  // namespace stan

#endif
