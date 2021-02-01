#ifndef STAN_MATH_OPENCL_REV_ADD_DIAG_HPP
#define STAN_MATH_OPENCL_REV_ADD_DIAG_HPP
#ifdef STAN_OPENCL
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/prim/add_diag.hpp>
#include <stan/math/opencl/prim/sum.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/adjoint_of.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>

namespace stan {
namespace math {

/**
 * Returns a Matrix with values added along the main diagonal
 *
 * @tparam T_m type of input kernel generator expression for the input matrix
 * @tparam T_a type of input kernel generator expression to add along the
 * diagonal
 *
 * @param mat input kernel generator expression
 * @param to_add scalar value or input kernel generator expression to add along
 * the diagonal
 * @return a kernel generator expressio with to_add added along main diagonal
 */
template <typename T_m, typename T_a,
          require_all_nonscalar_prim_or_rev_kernel_expression_t<T_m>* = nullptr,
          require_all_prim_or_rev_kernel_expression_t<T_a>* = nullptr,
          require_any_var_t<T_m, T_a>* = nullptr>
inline auto add_diag(const T_m& mat, const T_a& to_add) {
  const arena_t<T_m>& mat_arena = mat;
  const arena_t<T_a>& to_add_arena = to_add;

  var_value<matrix_cl<double>> res = add_diag(value_of(mat), value_of(to_add));

  reverse_pass_callback([mat_arena, to_add_arena, res]() mutable {
    if (!is_constant<T_m>::value) {
      adjoint_of(mat_arena) += res.adj();
    }
    if (!is_constant<T_a>::value) {
      if (!is_stan_scalar<T_a>::value) {
        auto& to_add_adj
            = forward_as<var_value<matrix_cl<double>>>(to_add_arena).adj();
        to_add_adj += diagonal(res.adj());
      } else {
        auto& to_add_adj = forward_as<var_value<double>>(to_add_arena).adj();
        to_add_adj += to_add_adj + sum(diagonal(res.adj()));
      }
    }
  });
  return res;
}
}  // namespace math
}  // namespace stan

#endif
#endif
