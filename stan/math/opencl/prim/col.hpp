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
 * This is equivalent to calling <code>m.col(i - 1)</code> and
 * assigning the resulting expression to a matrix_cl representing
 * a column vector.
 *
 * @tparam T type of the matrix
 * @param m input kernel generator expression.
 * @param j Column index (count from 1).
 * @return Specified column of the matrix.
 * @throw std::invalid_argument if j is out of range.
 */
template <typename T_a,
          typename = require_all_valid_expressions_and_none_scalar_t<T_a>>
inline auto col(T_a&& a, size_t j) {  // NOLINT
  return block(std::forward<T_a>(a), 0, j - 1, a.rows(), 1);
}
}  // namespace math
}  // namespace stan
#endif
#endif
