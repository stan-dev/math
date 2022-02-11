#ifndef STAN_MATH_OPENCL_PRIM_MULTIPLY_LOWER_TRI_SELF_TRANSPOSE_HPP
#define STAN_MATH_OPENCL_PRIM_MULTIPLY_LOWER_TRI_SELF_TRANSPOSE_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/prim/multiply.hpp>

namespace stan {
namespace math {

/**
 * Returns the result of multiplying the lower triangular
 * portion of the input matrix by its own transpose.
 *
 * @param x Matrix to multiply.
 * @return The lower triangular values in x times their own
 * transpose.
 */
template <typename T_x,
          require_all_kernel_expressions_and_none_scalar_t<T_x>* = nullptr>
inline auto multiply_lower_tri_self_transpose(T_x&& x) {
  matrix_cl<double> x_eval = std::forward<T_x>(x);
  if (x_eval.size() == 0) {
    return matrix_cl<double>(0, 0);
  }
  matrix_cl<double> x_lower(x_eval.buffer(), x_eval.rows(), x_eval.cols(),
                            matrix_cl_view::Lower);
  for (auto e : x_eval.write_events()) {
    x_lower.add_write_event(e);
  }
  matrix_cl<double> res = x_lower * transpose(x_lower);
  x_eval.add_read_event(x_lower.read_events().back());
  return res;
}
}  // namespace math
}  // namespace stan
#endif
#endif
