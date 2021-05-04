#ifndef STAN_MATH_PRIM_ERR_IS_LOWER_TRIANGULAR_HPP
#define STAN_MATH_PRIM_ERR_IS_LOWER_TRIANGULAR_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/err/is_not_nan.hpp>
#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {

/**
 * Return <code>true</code> is matrix is lower triangular.
 * A matrix x is not lower triangular if there is a non-zero entry
 * x[m, n] with m &lt; n. This function only inspect the upper and
 * triangular portion of the matrix, not including the diagonal.
 * @tparam EigMat A type derived from `EigenBase` with dynamic rows and columns
 * @param y Matrix to test
 * @return <code>true</code> is matrix is lower triangular
 */
template <typename EigMat, require_eigen_matrix_dynamic_t<EigMat>* = nullptr>
inline bool is_lower_triangular(const EigMat& y) {
  return y.unaryExpr([](auto&& x) { return is_not_nan(x) ? x : 1.0; })
      .transpose()
      .isUpperTriangular();
}

}  // namespace math
}  // namespace stan
#endif
