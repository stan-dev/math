#ifndef STAN_MATH_PRIM_FUN_LOG_DETERMINANT_LDLT_HPP
#define STAN_MATH_PRIM_FUN_LOG_DETERMINANT_LDLT_HPP

#include <stan/math/prim/fun/LDLT_factor.hpp>

namespace stan {
namespace math {

/**
 * Returns log(abs(det(A))) given a LDLT_factor of A
 *
 * @tparam T type of elements in the LDLT_factor
 * @param A LDLT_factor
 * @return the log(abs(det(A))
 */
template <typename T, require_not_rev_matrix_t<T>* = nullptr>
inline value_type_t<T> log_determinant_ldlt(LDLT_factor<T>& A) {
  if (A.matrix().size() == 0) {
    return 0;
  }

  return sum(log(A.ldlt().vectorD().array()));
}

}  // namespace math
}  // namespace stan

#endif
