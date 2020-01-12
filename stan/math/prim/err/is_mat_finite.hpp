#ifndef STAN_MATH_PRIM_ERR_IS_MAT_FINITE_HPP
#define STAN_MATH_PRIM_ERR_IS_MAT_FINITE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Return <code>true</code> is the specified matrix is finite.
 * @tparam T Scalar type of the matrix, requires class method
 * <code>.allFinite()</code>
 * @tparam R number of rows or Eigen::Dynamic
 * @tparam C number of columns or Eigen::Dynamic
 * @param y Matrix to test
 * @return <code>true</code> if the matrix is finite
 **/
template <typename T, int R, int C>
inline bool is_mat_finite(const Eigen::Matrix<T, R, C>& y) {
  return y.allFinite();
}

}  // namespace math
}  // namespace stan
#endif
