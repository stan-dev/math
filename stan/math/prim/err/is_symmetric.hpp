#ifndef STAN_MATH_PRIM_ERR_IS_SYMMETRIC_HPP
#define STAN_MATH_PRIM_ERR_IS_SYMMETRIC_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/constraint_tolerance.hpp>
#include <stan/math/prim/err/is_square.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return <code>true</code> if the matrix is square, and no element
 * not on the main diagonal is <code>NaN</code>.
 * @tparam EigMat A type derived from `EigenBase` with dynamic rows and columns
 * @param y Matrix to test
 * @return <code>true</code> if the matrix is square, and no
 *    element not on the main diagonal is <code>NaN</code>
 */
template <typename EigMat, require_eigen_matrix_dynamic_t<EigMat>* = nullptr>
inline bool is_symmetric(const EigMat& y) {
  const auto& y_ref = to_ref(y);
  if (!is_square(y_ref)) {
    return false;
  }
  Eigen::Index k = y_ref.rows();
  if (k == 1) {
    return true;
  }
  for (Eigen::Index m = 0; m < k; ++m) {
    for (Eigen::Index n = m + 1; n < k; ++n) {
      if (!(fabs(value_of(y_ref(m, n)) - value_of(y_ref(n, m)))
            <= CONSTRAINT_TOLERANCE)) {
        return false;
      }
    }
  }
  return true;
}

}  // namespace math
}  // namespace stan
#endif
