#ifndef STAN_MATH_OPENCL_PRIM_COL_HPP
#define STAN_MATH_OPENCL_PRIM_COL_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/kernel_generator.hpp>

namespace stan {
namespace math {

/** \ingroup opencl
 * Return the specified column of the specified kernel generator
 * expression using start-at-1 indexing.
 *
 * @tparam T_x type of input kernel generator expression a
 * @param x input kernel generator expression.
 * @param j Column index (count from 1).
 * @return Specified column of the matrix.
 * @throw std::invalid_argument if j is out of range.
 */
template <typename T_x,
          typename = require_all_kernel_expressions_and_none_scalar_t<T_x>>
inline auto col(T_x&& x, size_t j) {  // NOLINT
  return block(std::forward<T_x>(x), 0, j - 1, x.rows(), 1);
}
}  // namespace math
}  // namespace stan
#endif
#endif
