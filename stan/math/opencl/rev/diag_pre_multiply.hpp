#ifndef STAN_MATH_OPENCL_REV_DIAG_PRE_MULTIPLY_HPP
#define STAN_MATH_OPENCL_REV_DIAG_PRE_MULTIPLY_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/adjoint_of.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/eval.hpp>

namespace stan {
namespace math {

/**
 * Return the product of the diagonal matrix formed from the vector
 * or row_vector and a matrix.
 *
 * @tparam T1 type of the vector/row_vector
 * @tparam T2 type of the matrix
 * @param v1 input vector/row_vector
 * @param v2 input matrix
 *
 * @return product of the diagonal matrix formed from the
 * vector or row_vector and a matrix.
 */
template <
    typename T1, typename T2, require_any_var_t<T1, T2>* = nullptr,
    require_all_nonscalar_prim_or_rev_kernel_expression_t<T1, T2>* = nullptr>
inline var_value<matrix_cl<double>> diag_pre_multiply(T1&& v1, T2&& v2) {
  arena_t<T1> v1_arena = std::forward<T1>(v1);
  arena_t<T2> v2_arena = std::forward<T2>(v2);

  matrix_cl<double> res_val
      = diag_pre_multiply(value_of(v1_arena), value_of(v2_arena));

  return make_callback_var(
      res_val,
      [v1_arena, v2_arena](const vari_value<matrix_cl<double>>& res) mutable {
        if (v1_arena.cols() == 1) {
          if (!is_constant<std::decay_t<T1>>::value) {
            adjoint_of(v1_arena)
                += rowwise_sum(elt_multiply(res.adj(), value_of(v2_arena)));
          }
          if (!is_constant<std::decay_t<T2>>::value) {
            adjoint_of(v2_arena) += elt_multiply(
                res.adj(), rowwise_broadcast(value_of(v1_arena)));
          }
        } else {
          if (!is_constant<std::decay_t<T1>>::value) {
            adjoint_of(transpose(v1_arena))
                += rowwise_sum(elt_multiply(res.adj(), value_of(v2_arena)));
          }
          if (!is_constant<std::decay_t<T2>>::value) {
            adjoint_of(v2_arena) += elt_multiply(
                res.adj(), rowwise_broadcast(transpose(value_of(v1_arena))));
          }
        }
      });
}

}  // namespace math
}  // namespace stan

#endif
#endif
