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
template <typename T, bool alloc_in_arena,
          require_not_rev_matrix_t<T>* = nullptr>
inline value_type_t<T> log_determinant_ldlt(
    const LDLT_factor<T, alloc_in_arena>& A) {
  if (A.matrix().size() == 0) {
    return 0;
  }

  return sum(log(A.ldlt().vectorD().array()));
}

}  // namespace math
}  // namespace stan

#endif
