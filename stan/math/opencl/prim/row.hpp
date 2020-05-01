#ifndef STAN_MATH_OPENCL_PRIM_ROW_HPP
#define STAN_MATH_OPENCL_PRIM_ROW_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/kernel_generator.hpp>

namespace stan {
namespace math {

/** \ingroup opencl
 * Return the specified row of the specified kernel generator
 * expression using start-at-1 indexing.
 *
 * This is equivalent to calling <code>m.row(i - 1)</code> and
 * assigning the resulting expression to a matrix_cl representing
 * a row vector.
 *
 * @tparam T type of the matrix
 * @param a input kernel generator expression.
 * @param j Row index (count from 1).
 * @return Specified row of the matrix.
 * @throw std::out_of_range if j is out of range.
 */
template <typename T_a,
          typename = require_all_valid_expressions_and_none_scalar_t<T_a>>
inline auto row(T_a&& a, size_t j) {  // NOLINT
  return block(std::forward<T_a>(a), j - 1, 0, 1, a.cols());
}
}  // namespace math
}  // namespace stan
#endif
#endif
