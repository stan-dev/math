#ifndef STAN_MATH_REV_FUN_FILL_HPP
#define STAN_MATH_REV_FUN_FILL_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Fill the specified container with the specified value.
 *
 * The specified matrix is filled by element.
 *
 * @tparam VarMat a `var_value` with inner type from `EigenBase`
 * @tparam Var A var.
 *
 * @param x Container.
 * @param y Value.
 */
template <typename VarMat, typename Var, require_var_matrix_t<VarMat>* = nullptr,
          require_var_t<Var>* = nullptr>
inline void fill(VarMat& x, const Var& y) {
  arena_t<plain_type_t<value_type_t<VarMat>>> prev_vals = x.val().eval();
  x.vi_->val_.fill(y.val());
  reverse_pass_callback([x, y, prev_vals]() mutable {
    x.vi_->val_ = prev_vals;
    y.adj() += x.adj().sum();
    x.adj().setZero();
  });
}

/**
 * Fill the specified container with the specified value.
 *
 * The specified matrix is filled by element.
 *
 * @tparam VarMat a `var_value` with inner type from `EigenBase`
 * @tparam Arith An arithmetic type.
 *
 * @param x Container.
 * @param y Value.
 */
template <typename VarMat, typename Arith, require_var_matrix_t<VarMat>* = nullptr,
          require_arithmetic_t<Arith>* = nullptr>
inline void fill(VarMat& x, const Arith& y) {
  arena_t<plain_type_t<value_type_t<VarMat>>> prev_vals = x.val().eval();
  x.vi_->val_.fill(y);
  reverse_pass_callback([x, prev_vals]() mutable {
    x.vi_->val_ = prev_vals;
    x.adj().setZero();
  });
}

}  // namespace math
}  // namespace stan

#endif
