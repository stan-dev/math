#ifndef STAN_MATH_OPENCL_REV_ROWS_DOT_PRODUCT_HPP
#define STAN_MATH_OPENCL_REV_ROWS_DOT_PRODUCT_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/prim/rows_dot_product.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/adjoint_of.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/prim/fun/size_zero.hpp>

namespace stan {
namespace math {

/**
 * Returns the dot product of rows of the specified matrices.
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
inline var_value<matrix_cl<double>> rows_dot_product(T1&& v1, T2&& v2) {
  check_matching_sizes("rows_dot_product(OpenCL)", "v1", v1, "v2", v2);

  arena_t<T1> v1_arena = std::forward<T1>(v1);
  arena_t<T2> v2_arena = std::forward<T2>(v2);

  matrix_cl<double> res_val = rows_dot_product(value_of(v1), value_of(v2));

  return make_callback_var(
      res_val,
      [v1_arena, v2_arena](const vari_value<matrix_cl<double>>& res) mutable {
        auto v1_deriv
            = elt_multiply(rowwise_broadcast(res.adj()), value_of(v2_arena));
        auto v2_deriv
            = elt_multiply(rowwise_broadcast(res.adj()), value_of(v1_arena));
        results(adjoint_of(v1_arena), adjoint_of(v2_arena))
            += expressions(calc_if<is_var<T1>::value>(v1_deriv),
                           calc_if<is_var<T2>::value>(v2_deriv));
      });
}

}  // namespace math
}  // namespace stan

#endif
#endif
