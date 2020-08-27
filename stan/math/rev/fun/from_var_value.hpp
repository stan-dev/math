#ifndef STAN_MATH_REV_FUN_FROM_VAR_VALUE_HPP
#define STAN_MATH_REV_FUN_FROM_VAR_VALUE_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/functor/arena_matrix.hpp>
#include <stan/math/rev/functor/reverse_pass_callback.hpp>

namespace stan {
namespace math {

/**
 * Converts `var_value` into an Eigen Matrix. Adjoint is propagated back to
 * arguument in the reverse pass.
 *
 * @tparam T type of the input
 * @param a matrix to convert
 */
template <typename T, require_var_vt<is_eigen, T>* = nullptr>
Eigen::Matrix<var, value_type_t<T>::RowsAtCompileTime,
              value_type_t<T>::ColsAtCompileTime>
from_var_value(const T& a) {
  arena_matrix<Eigen::Matrix<var, value_type_t<T>::RowsAtCompileTime,
                             value_type_t<T>::ColsAtCompileTime>>
      res(a.val());
  reverse_pass_callback([res, a]() mutable { a.vi_->adj_ += res.adj(); });
  return res;
}

}  // namespace math
}  // namespace stan

#endif
