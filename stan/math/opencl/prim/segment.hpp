#ifndef STAN_MATH_OPENCL_PRIM_SEGMENT_HPP
#define STAN_MATH_OPENCL_PRIM_SEGMENT_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/err.hpp>
#include <stan/math/opencl/prim/block.hpp>

namespace stan {
namespace math {

/**
 * Return the specified number of elements as a row/column vector starting
 * from the specified element - 1 of the specified row/column vector.
 *
 * @tparam T_x type of input kernel generator expression x
 * @param x input kernel generator expression.
 * @param i Starting row/column + 1.
 * @param n Number of rows/columns in segment.
 * @throw std::out_of_range if either index is out of range.
 */
template <typename T_x,
          require_nonscalar_prim_or_rev_kernel_expression_t<T_x>* = nullptr>
inline auto segment(T_x&& x, size_t i, size_t n) {
  check_vector("segment (OpenCL)", "x", x);
  if (x.rows() == 1) {
    return block(std::forward<T_x>(x), 1, i, 1, n);
  } else {
    return block(std::forward<T_x>(x), i, 1, n, 1);
  }
}
}  // namespace math
}  // namespace stan
#endif
#endif
