#ifndef STAN_MATH_OPENCL_PRIM_SUB_ROW_HPP
#define STAN_MATH_OPENCL_PRIM_SUB_ROW_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/prim/block.hpp>

namespace stan {
namespace math {

/**
 * Return a 1 x ncols subrow starting at (i-1, j-1).
 *
 * @tparam T_x type of input kernel generator expression x
 * @param x input kernel generator expression.
 * @param i Starting row + 1.
 * @param j Starting column + 1.
 * @param ncols Number of columns in block.
 * @throw std::out_of_range if either index is out of range.
 */
template <typename T_x,
          typename = require_nonscalar_prim_or_rev_kernel_expression_t<T_x>>
inline auto sub_row(T_x&& x, size_t i, size_t j, size_t ncols) {
  return block(std::forward<T_x>(x), i, j, 1, ncols);
}
}  // namespace math
}  // namespace stan
#endif
#endif
