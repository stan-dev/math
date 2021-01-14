#ifndef STAN_MATH_OPENCL_REV_COLUMNS_DOT_SELF_HPP
#define STAN_MATH_OPENCL_REV_COLUMNS_DOT_SELF_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/adjoint_of.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/prim/fun/size_zero.hpp>

namespace stan {
namespace math {

/**
 * Returns the dot product of each column of a matrix with itself.
 *
 * @tparam T type of the matrix
 * @param v Matrix.
 */
template <typename T,
          require_var_vt<is_kernel_expression_and_not_scalar, T>* = nullptr>
inline var_value<matrix_cl<double>> columns_dot_self(T&& v) {
  if (size_zero(v)) {
    return var_value<matrix_cl<double>>(constant(0.0, 1, v.cols()));
  }

  arena_t<T> v_arena;
  if ((std::is_rvalue_reference<T&&>::value && is_matrix_cl<T>::value)
      || is_var<T>::value) {
    v_arena = std::forward<T>(v);
  }

  matrix_cl<double> res_val;
  results(v_arena, res_val) = expressions(
      calc_if<((std::is_lvalue_reference<T>::value || !is_matrix_cl<T>::value)
               && !is_var<T>::value)>(value_of(v)),
      colwise_sum(square(value_of(v))));
  while (res_val.rows() > 1) {
    res_val = colwise_sum(res_val).eval();
  }

  return make_callback_var(
      res_val, [v_arena](const vari_value<matrix_cl<double>>& res) mutable {
        v_arena.adj() = v_arena.adj()
                        + elt_multiply(colwise_broadcast(res.adj() * 2),
                                       value_of(v_arena));
      });
}

}  // namespace math
}  // namespace stan

#endif
#endif
