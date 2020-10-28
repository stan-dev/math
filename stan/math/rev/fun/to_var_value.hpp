#ifndef STAN_MATH_REV_FUN_TO_VAR_VALUE_HPP
#define STAN_MATH_REV_FUN_TO_VAR_VALUE_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/core/arena_matrix.hpp>
#include <stan/math/rev/functor/reverse_pass_callback.hpp>

namespace stan {
namespace math {

/**
 * Converts an Eigen matrix (or vector or row_vector) or expression of `var`s
 * into `var_value`. Adjoint is propagated back to argument in the reverse
 * pass.
 *
 * @tparam T type of the input
 * @param a matrix to convert
 */
template <typename T, require_eigen_vt<is_var, T>* = nullptr>
var_value<Eigen::Matrix<double, T::RowsAtCompileTime, T::ColsAtCompileTime>>
to_var_value(const T& a) {
  arena_matrix<plain_type_t<T>> a_arena = a;
  var_value<promote_scalar_t<double, T>> res(a_arena.val());
  reverse_pass_callback(
      [res, a_arena]() mutable { a_arena.adj() += res.adj(); });
  return res;
}

}  // namespace math
}  // namespace stan

#endif
