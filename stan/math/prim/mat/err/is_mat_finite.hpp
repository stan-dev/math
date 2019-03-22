#ifndef STAN_MATH_PRIM_MAT_ERR_IS_MAT_FINITE_HPP
#define STAN_MATH_PRIM_MAT_ERR_IS_MAT_FINITE_HPP

#include <stan/math/prim/mat/fun/value_of.hpp>
#include <Eigen/Dense>
#include <boost/math/special_functions/fpclassify.hpp>

namespace stan {
namespace math {

/**
 * Return <code>true</code> is the specified matrix is finite.
 * @tparam T Scalar type of the matrix, requires class method
 *   <code>.size()</code>
 * @tparam R Compile time rows of the matrix
 * @tparam C Compile time columns of the matrix
 * @param y Matrix to test
 * @return <code>true</code> if the matrix is finite
 **/
template <typename T, int R, int C>
inline bool is_mat_finite(const Eigen::Matrix<T, R, C>& y) {
  if (!value_of(y).allFinite()) {
    for (int n = 0; n < y.size(); ++n) {
      if (!(boost::math::isfinite)(y(n)))
        return false;
    }
  }
  return true;
}

}  // namespace math
}  // namespace stan
#endif
