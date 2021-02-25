#ifndef STAN_MATH_OPENCL_PRIM_REVERSE_HPP
#define STAN_MATH_OPENCL_PRIM_REVERSE_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/kernel_generator.hpp>

namespace stan {
namespace math {

/**
 * Return reversed view into the specified vector or row vector.
 *
 * @tparam T type of expression
 * @param x vector or row vector to reverse
 * @return reversed vector or row vector
 */
template <typename T_x,
          require_all_kernel_expressions_and_none_scalar_t<T_x>* = nullptr>
inline auto reverse(T_x&& x) {
  check_vector("reverse(OpenCL)", "x", x);
  return indexing(std::forward<T_x>(x),
                  x.rows() - row_index(x.rows(), x.cols()) - 1,
                  x.cols() - col_index() - 1);
}
}  // namespace math
}  // namespace stan
#endif
#endif
