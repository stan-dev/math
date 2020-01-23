#ifndef STAN_MATH_PRIM_FUN_SD_HPP
#define STAN_MATH_PRIM_FUN_SD_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/variance.hpp>
#include <stan/math/prim/fun/sqrt.hpp>
#include <vector>
#include <cmath>

namespace stan {
namespace math {

/**
 * Returns the unbiased sample standard deviation of the
 * coefficients in the specified column vector.
 *
 * @tparam T type of elements in the vector
 * @param v Specified vector.
 * @return Sample variance of vector.
 */
template <typename T>
inline return_type_t<T> sd(const std::vector<T>& v) {
  check_nonzero_size("sd", "v", v);
  if (v.size() == 1) {
    return 0.0;
  }
  return sqrt(variance(v));
}

/**
 * Returns the unbiased sample standard deviation of the
 * coefficients in the specified vector, row vector, or matrix.
 *
 * @tparam T type of elements in the vector, row vector, or matrix
 * @tparam R number of rows, can be Eigen::Dynamic
 * @tparam C number of columns, can be Eigen::Dynamic
 *
 * @param m Specified vector, row vector or matrix.
 * @return Sample variance.
 */
template <typename T, int R, int C>
inline return_type_t<T> sd(const Eigen::Matrix<T, R, C>& m) {
  using std::sqrt;
  check_nonzero_size("sd", "m", m);
  if (m.size() == 1) {
    return 0.0;
  }
  return sqrt(variance(m));
}

}  // namespace math
}  // namespace stan

#endif
