#ifndef STAN_MATH_OPENCL_REV_DIAG_POST_MULTIPLY_HPP
#define STAN_MATH_OPENCL_REV_DIAG_POST_MULTIPLY_HPP
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
 * Return the product of a matrix and the diagonal matrix formed from the vector
 * or row_vector.
 *
 * @tparam T1 type of the matrix
 * @tparam T2 type of the vector/row_vector
 * @param v1 input matrix
 * @param v2 input vector/row_vector
 *
 * @return product of a matrix and the diagonal matrix formed from the
 * vector or row_vector.
 */
template <
    typename T1, typename T2, require_any_var_t<T1, T2>* = nullptr,
    require_all_nonscalar_prim_or_rev_kernel_expression_t<T1, T2>* = nullptr>
inline var_value<matrix_cl<double>> diag_post_multiply(T1&& v1, T2&& v2) {
  arena_t<T1> v1_arena = std::forward<T1>(v1);
  arena_t<T2> v2_arena = std::forward<T2>(v2);

  matrix_cl<double> res_val
      = diag_post_multiply(value_of(v1_arena), value_of(v2_arena));

  return make_callback_var(
      res_val,
      [v1_arena, v2_arena](const vari_value<matrix_cl<double>>& res) mutable {
        if (v2_arena.rows() == 1) {
          auto v1_adj_inc
              = elt_multiply(res.adj(), colwise_broadcast(value_of(v2_arena)));
          auto v2_adj_inc
              = colwise_sum(elt_multiply(res.adj(), value_of(v1_arena)));
          matrix_cl<double> tmp;
          auto&& v1_adj = adjoint_of(v1_arena);
          results(v1_adj, tmp) = expressions(
              calc_if<!is_constant<std::decay_t<T1>>::value>(v1_adj
                                                             + v1_adj_inc),
              calc_if<!is_constant<std::decay_t<T2>>::value>(v2_adj_inc));

          if (!is_constant<std::decay_t<T2>>::value) {
            while (tmp.rows() > 1) {
              tmp = eval(colwise_sum(tmp));
            }
            adjoint_of(v2_arena) += tmp;
          }
        } else {
          auto v1_adj_inc = elt_multiply(
              res.adj(), colwise_broadcast(transpose(value_of(v2_arena))));
          auto v2_adj_inc
              = colwise_sum(elt_multiply(res.adj(), value_of(v1_arena)));
          matrix_cl<double> tmp;
          auto&& v1_adj = adjoint_of(v1_arena);
          results(v1_adj, tmp) = expressions(
              calc_if<!is_constant<std::decay_t<T1>>::value>(v1_adj
                                                             + v1_adj_inc),
              calc_if<!is_constant<std::decay_t<T2>>::value>(v2_adj_inc));
          if (!is_constant<std::decay_t<T2>>::value) {
            while (tmp.rows() > 1) {
              tmp = eval(colwise_sum(tmp));
            }
            adjoint_of(v2_arena) += transpose(tmp);
          }
        }
      });
}

}  // namespace math
}  // namespace stan

#endif
#endif
