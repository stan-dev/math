#ifndef DIAG_PRE_MULTIPLY_HPP
#define DIAG_PRE_MULTIPLY_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/adjoint_of.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/prim/fun/size_zero.hpp>

namespace stan {
namespace math {

/**
 * Returns the dot product of columns of the specified matrices.
 *
 * @tparam T1 type of the first matrix
 * @tparam T2 type of the second matrix
 *
 * @param v1 Matrix of first vectors.
 * @param v2 Matrix of second vectors.
 * @return Dot product of the vectors.
 * @throw std::invalid_argument If the vectors are not the same
 * size
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
        if (!is_constant<std::decay_t<T1>>::value) {
          auto v1_adj
              = forward_as<var_value<matrix_cl<double>>>(v1_arena).adj();
          v1_adj = v1_adj
                   + rowwise_sum(elt_multiply(res.adj(), value_of(v2_arena)));
        }
        if (!is_constant<std::decay_t<T2>>::value) {
          auto v2_adj
              = forward_as<var_value<matrix_cl<double>>>(v2_arena).adj();
          v2_adj = v2_adj
                   + elt_multiply(res.adj(),
                                  rowwise_broadcast(value_of(v1_arena)));
        }
      });
}

}  // namespace math
}  // namespace stan

#endif
#endif
