#ifndef STAN_MATH_REV_META_RETURN_VAR_MATRIX
#define STAN_MATH_REV_META_RETURN_VAR_MATRIX

#include <stan/math/rev/meta/is_var.hpp>
#include <stan/math/prim/meta.hpp>

namespace stan {

/**
 * Given an Eigen type and several inputs, determine if a matrix should be
 * `var<Matrix>` or `Matrix<var>`.
 *
 * `Matrix` will be a plain type
 *
 * @tparam ReturnType An Eigen matrix used for composing the `var<Matrix>` or
 *  `Matrix<var>` type.
 * @tparam Types Parameter pack holding any mix of types. If any of `Types`
 *  are a `var<Matrix>` this holds a `var<Matrix>` type.
 *  Else the type will be `Matrix<var>`
 */
template <typename ReturnType, typename... Types>
using return_var_matrix_t = std::conditional_t<
    is_any_var_matrix<ReturnType, Types...>::value,
    stan::math::var_value<
        stan::math::promote_scalar_t<double, plain_type_t<ReturnType>>>,
    stan::math::promote_scalar_t<stan::math::var_value<double>,
                                 plain_type_t<ReturnType>>>;
}  // namespace stan

#endif
