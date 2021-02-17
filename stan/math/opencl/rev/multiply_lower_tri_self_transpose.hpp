#ifndef STAN_MATH_OPENCL_REV_MULTIPLY_LOWER_TRI_SELF_TRANSPOSE_HPP
#define STAN_MATH_OPENCL_REV_MULTIPLY_LOWER_TRI_SELF_TRANSPOSE_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/prim/multiply_lower_tri_self_transpose.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>

namespace stan {
namespace math {

/**
 * Returns the result of multiplying the lower triangular
 * portion of the input matrix by its own transpose.
 *
 * @param A Matrix to multiply.
 * @return The lower triangular values in L times their own
 * transpose.
 */
template <typename T,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr>
inline var_value<matrix_cl<double>> multiply_lower_tri_self_transpose(
    const var_value<T>& A) {
  return make_callback_var(
      multiply_lower_tri_self_transpose(A.val()),
      [A](vari_value<matrix_cl<double>>& res) mutable {
        if (A.size() != 0) {
          matrix_cl<double> A_val_lower(A.val().buffer(), A.val().rows(),
                                        A.val().cols(), matrix_cl_view::Lower);
          for (auto e : A.val().write_events()) {
            A_val_lower.add_write_event(e);
          }
          matrix_cl<double> prod
              = (res.adj() + transpose(res.adj())) * A_val_lower;
          prod.view(matrix_cl_view::Lower);
          A.adj() += prod;
        }
      });
}

}  // namespace math
}  // namespace stan

#endif
#endif
