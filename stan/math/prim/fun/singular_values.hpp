#ifndef STAN_MATH_PRIM_FUN_SINGULAR_VALUES_HPP
#define STAN_MATH_PRIM_FUN_SINGULAR_VALUES_HPP

#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Return the vector of the singular values of the specified matrix
 * in decreasing order of magnitude.
 * <p>See the documentation for <code>svd()</code> for
 * information on the singular values.
 *
 * @tparam T type of elements in the matrix
 * @param m Specified matrix.
 * @return Singular values of the matrix.
 */
template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> singular_values(
    const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& m) {
  if (m.size() == 0) {
    return {};
  }

  return Eigen::JacobiSVD<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >(m)
      .singularValues();
}

}  // namespace math
}  // namespace stan

#endif
