#ifndef STAN_MATH_OPENCL_PRIM_TAIL_HPP
#define STAN_MATH_OPENCL_PRIM_TAIL_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/prim/block.hpp>
#include <stan/math/prim/err.hpp>

namespace stan {
namespace math {

/**
 * Return the specified number of elements as a vector or row vector (same as
 * input) from the back of the specified vector or row vector.
 *
 * @tparam T_x type of input kernel generator expression.
 * @param x input kernel generator expression.
 * @param n Size of return.
 * @return The first n elements of v.
 * @throw std::out_of_range if n is out of range.
 */
template <typename T_x,
          require_nonscalar_prim_or_rev_kernel_expression_t<T_x>* = nullptr>
inline auto tail(T_x&& x, size_t n) {  // NOLINT
  check_vector("head (OpenCL)", "x", x);
  if (n != 0) {
    check_vector_index("head", "n", x, n);
  }
  if (x.rows() == 1) {
    return block(std::forward<T_x>(x), 1, 1 + x.cols() - n, 1, n);
  } else {
    return block(std::forward<T_x>(x), 1 + x.rows() - n, 1, n, 1);
  }
}
}  // namespace math
}  // namespace stan
#endif
#endif
