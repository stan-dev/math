#ifndef STAN_MATH_REV_FUN_IDENTITY_CONSTRAIN_HPP
#define STAN_MATH_REV_FUN_IDENTITY_CONSTRAIN_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/identity_constrain.hpp>
#include <stan/math/rev/fun/to_var_value.hpp>
namespace stan {
namespace math {

/**
 * Returns the result of applying the identity constraint
 * transform to the input. This specialization handles convert Eigen matrices
 * of doubles to var matrix types.
 *
 * @tparam ArithEigMat Any type.
 * @tparam Types Any type with one of `ArithEigMat` and `Types` being a `var_value`
 * matrix.
 * @param[in] x object
 * @return transformed input
 */
template <typename ArithEigMat, typename... Types,
          require_eigen_vt<std::is_arithmetic, ArithEigMat>* = nullptr,
          require_any_var_matrix_t<ArithEigMat, Types...>* = nullptr>
inline auto identity_constrain(ArithEigMat&& x, Types&&... /* args */) {
  return var_value<plain_type_t<ArithEigMat>>(x);
}

template <typename VarEigMat, typename... Types, require_eigen_vt<is_var, VarEigMat>* = nullptr,
          require_any_var_matrix_t<VarEigMat, Types...>* = nullptr>
inline auto identity_constrain(VarEigMat&& x, Types&&... /* args */) {
  return to_var_value(x);
}

template <typename VarMat, typename... Types, require_var_matrix_t<VarMat>* = nullptr,
          require_any_var_matrix_t<VarMat, Types...>* = nullptr>
inline auto identity_constrain(VarMat&& x, Types&&... /* args */) {
  return x;
}

}  // namespace math
}  // namespace stan

#endif
