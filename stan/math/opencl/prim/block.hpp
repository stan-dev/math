#ifndef STAN_MATH_OPENCL_PRIM_BLOCK_HPP
#define STAN_MATH_OPENCL_PRIM_BLOCK_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/kernel_generator.hpp>

namespace stan {
namespace math {

/**
 * Return a nrows x ncols submatrix starting at (i-1, j-1).
 *
 * @tparam T type of elements in the matrix
 * @param x Matrix.
 * @param i Starting row.
 * @param j Starting column.
 * @param nrows Number of rows in block.
 * @param ncols Number of columns in block.
 * @throw std::out_of_range if either index is out of range.
 */
template <typename T_x,
          typename = require_all_kernel_expressions_and_none_scalar_t<T_x>>
inline auto block(T_x&& x, size_t i, size_t j, size_t nrows,
                  size_t ncols) {  // NOLINT
  return block_zero_based(std::forward<T_x>(x), i - 1, j - 1, nrows, ncols);
}
}  // namespace math
}  // namespace stan
#endif
#endif
