#ifndef STAN_MATH_OPENCL_PRIM_ROW_HPP
#define STAN_MATH_OPENCL_PRIM_ROW_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/prim/block.hpp>

namespace stan {
namespace math {

/** \ingroup opencl
 * Return the specified row of the specified kernel generator
 * expression using start-at-1 indexing.
 *
 * @tparam T_x type of input kernel generator expression x
 * @param x input kernel generator expression.
 * @param j Row index (count from 1).
 * @return Specified row of the matrix.
 * @throw std::out_of_range if j is out of range.
 */
template <typename T_x,
          typename = require_nonscalar_prim_or_rev_kernel_expression_t<T_x>>
inline auto row(T_x&& x, size_t j) {  // NOLINT
  return block(std::forward<T_x>(x), j, 1, 1, x.cols());
}
}  // namespace math
}  // namespace stan
#endif
#endif
