#ifndef STAN_MATH_PRIM_FUN_SINGULAR_VALUES_HPP
#define STAN_MATH_PRIM_FUN_SINGULAR_VALUES_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Return the vector of the singular values of the specified matrix
 * in decreasing order of magnitude.
 * <p>See the documentation for <code>svd()</code> for
 * information on the singular values.
 *
 * @tparam EigMat type of the matrix
 * @param m Specified matrix.
 * @return Singular values of the matrix.
 */
template <typename EigMat, require_eigen_matrix_t<EigMat>* = nullptr>
Eigen::Matrix<value_type_t<EigMat>, Eigen::Dynamic, 1> singular_values(
    const EigMat& m) {
  if (m.size() == 0) {
    return {};
  }

  return Eigen::JacobiSVD<Eigen::Matrix<value_type_t<EigMat>, Eigen::Dynamic,
                                        Eigen::Dynamic> >(m)
      .singularValues();
}

}  // namespace math
}  // namespace stan

#endif
