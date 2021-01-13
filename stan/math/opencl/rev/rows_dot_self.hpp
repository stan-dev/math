#ifndef STAN_MATH_OPENCL_REV_ROWS_DOT_SELF_HPP
#define STAN_MATH_OPENCL_REV_ROWS_DOT_SELF_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/prim/rows_dot_self.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/adjoint_of.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/prim/fun/size_zero.hpp>

namespace stan {
namespace math {

/**
 * Returns the dot product of each row of a matrix with itself.
 *
 * @tparam T type of the matrix
 * @param v Matrix.
 */
template <typename T,
          require_var_vt<is_kernel_expression_and_not_scalar, T>* = nullptr>
inline var_value<matrix_cl<double>> rows_dot_self(T&& v) {
  arena_t<T> v_arena = std::forward<T>(v);

  matrix_cl<double> res_val = rows_dot_self(value_of(v));

  return make_callback_var(
      res_val, [v_arena](const vari_value<matrix_cl<double>>& res) mutable {
        v_arena.adj() = v_arena.adj()
                        + elt_multiply(rowwise_broadcast(res.adj() * 2),
                                       value_of(v_arena));
      });
}

}  // namespace math
}  // namespace stan

#endif
#endif
