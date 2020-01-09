#ifndef STAN_MATH_PRIM_FUN_LOG_DETERMINANT_LDLT_HPP
#define STAN_MATH_PRIM_FUN_LOG_DETERMINANT_LDLT_HPP

#include <stan/math/prim/fun/LDLT_factor.hpp>

namespace stan {
namespace math {

/**
 * Returns log(abs(det(A))) given a LDLT_factor of A
 *
 * @tparam T type of elements in the LDLT_factor
 * @tparam R number of rows, can be Eigen::Dynamic
 * @tparam C number of columns, can be Eigen::Dynamic
 *
 * @param A LDLT_factor
 * @return the log(abs(det(A))
 */
template <int R, int C, typename T>
inline T log_determinant_ldlt(LDLT_factor<T, R, C> &A) {
  if (A.rows() == 0) {
    return 0;
  }

  return A.log_abs_det();
}

}  // namespace math
}  // namespace stan

#endif
