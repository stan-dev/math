#ifndef STAN_MATH_PRIM_MAT_ERR_IS_CHOLESKY_FACTOR_HPP
#define STAN_MATH_PRIM_MAT_ERR_IS_CHOLESKY_FACTOR_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/err/is_positive.hpp>
#include <stan/math/prim/mat/err/is_lower_triangular.hpp>
#include <stan/math/prim/scal/err/is_less_or_equal.hpp>

namespace stan {
namespace math {

/**
 * Return <code>true</code> if y is a valid Choleksy factor, if
 * number of rows is not less than the number of columns, if there
 * are no 0 columns, and no element in matrix is <code>NaN</code>.
 * A Cholesky factor is a lower triangular matrix whose diagonal
 * elements are all positive.  Note that Cholesky factors need not
 * be square, but requires at least as many rows M as columns N
 * (i.e., M &gt;= N).
 * @tparam T_y Type of elements of Cholesky factor, requires class method
 *   <code>.rows()</code>, <code>.cols()</code>, and <code>.diagonal()</code>
 * @param y Matrix to test
 * @return <code>true</code> if y is a valid Choleksy factor, if
 *   number of rows is not less than the number of columns,
 *   if there are no 0 columns, and no element in matrix is <code>NaN</code>
 */
template <typename T_y, enable_if_eigen<T_y>* = nullptr>
inline bool is_cholesky_factor(const T_y& y) {
  return is_less_or_equal(y.cols(), y.rows()) && is_positive(y.cols())
         && is_lower_triangular(y) && is_positive(y.diagonal());
}

}  // namespace math
}  // namespace stan
#endif
