#ifndef STAN_MATH_PRIM_FUN_PROMOTE_SCALAR_HPP
#define STAN_MATH_PRIM_FUN_PROMOTE_SCALAR_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/promote_scalar.hpp>

namespace stan {
namespace math {

/**
 * Struct to hold static function for promoting underlying scalar
 * types.  This specialization is for `var_value<Eigen::MatrixXd>` inputs.
 *
 * @tparam T return scalar type
 * @tparam S input matrix or vector or row vector type for static nested
 * function, which must have a scalar type assignable to T
 */
template <typename T, typename S>
struct promote_scalar_struct<T, S, require_var_matrix_t<S>> {
  /**
   * Return the matrix consisting of the recursive promotion of
   * the elements of the input matrix to the scalar type specified
   * by the return template parameter.
   *
   * @param x input matrix.
   * @return matrix with values promoted from input vector.
   */
  static auto apply(const S& x) { return x; }
};

/**
 * Struct to hold static function for promoting underlying scalar
 * types.  This specialization is for Eigen inputs.
 *
 * @tparam T return scalar type
 * @tparam S input matrix or vector or row vector type for static nested
 * function, which must have a scalar type assignable to T
 */
template <typename T, typename S>
struct promote_scalar_struct<T, S, require_all_t<is_var_matrix<T>, is_eigen<S>>>> {
  /**
   * Return the matrix consisting of the recursive promotion of
   * the elements of the input matrix to the scalar type specified
   * by the return template parameter.
   *
   * @param x input matrix.
   * @return matrix with values promoted from input vector.
   */
  static auto apply(const S& x) { return T(x); }
};


}  // namespace math
}  // namespace stan

#endif
