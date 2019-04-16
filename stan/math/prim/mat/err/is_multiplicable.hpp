#ifndef STAN_MATH_PRIM_MAT_ERR_IS_MULTIPLICABLE_HPP
#define STAN_MATH_PRIM_MAT_ERR_IS_MULTIPLICABLE_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/scal/err/is_size_match.hpp>
#include <stan/math/prim/scal/err/is_positive.hpp>

namespace stan {
namespace math {

/**
 * Return <code>true</code> if the matrices can be multiplied.
 * This checks the runtime sizes to determine whether the two
 * matrices are multiplicable. This allows Eigen matrices,
 * vectors, and row vectors to be checked.
 * @tparam T1 Type of first matrix, requires class method <code>.rows()</code>
 *   and <code>.cols()</code>
 * @tparam T2 Type of second matrix, requires class method <code>.rows()</code>
 *   and <code>.cols()</code>
 * @param y1 First matrix to test
 * @param y2 Second matrix to test
 * @return <code>true</code> if the matrices are not multiplicable and if
 *   either matrix is size 0 for either rows and columns
 **/
template <typename T1, typename T2>
inline bool is_multiplicable(const T1& y1, const T2& y2) {
  return (is_positive(y1.rows()) && is_positive(y2.cols())
          && is_size_match(y1.cols(), y2.rows()) && is_positive(y1.cols()));
}

}  // namespace math
}  // namespace stan
#endif
