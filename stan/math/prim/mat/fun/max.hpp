#ifndef STAN_MATH_PRIM_MAT_FUN_MAX_HPP
#define STAN_MATH_PRIM_MAT_FUN_MAX_HPP

#include <stan/math/prim/arr/err/check_nonzero_size.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <algorithm>
#include <limits>
#include <vector>

namespace stan {
namespace math {

/**
 * Returns the maximum coefficient in the specified
 * column vector.
 * @param x specified vector
 * @return maximum coefficient value in the vector
 * @throws <code>std::invalid_argument</code> if the vector is size zero
 */
 */
inline int max(const std::vector<int>& x) {
  check_nonzero_size("max", "int vector", x);
  Eigen::Map<const Eigen::Matrix<int, Eigen::Dynamic, 1>> m(&x[0], x.size());
  return m.maxCoeff();
}

/**
 * Returns the maximum coefficient in the specified
 * column vector.
 * @tparam type of values being compared and returned
 * @param x specified vector
 * @return maximum coefficient value in the vector, or -infinity if it is size zero
 */
template <typename T>
inline T max(const std::vector<T>& x) {
  if (x.size() == 0) {
    return -std::numeric_limits<T>::infinity();
  }
  Eigen::Map<const Eigen::Matrix<T, Eigen::Dynamic, 1>> m(&x[0], x.size());
  return m.maxCoeff();
}

/**
 * Returns the maximum coefficient in the specified
 * matrix, vector, or row vector.
 * @tparam T type of values being compared and returned
 * @tparam R number of rows or Eigen::Dynamic
 * @tparam C number of columns or Eigen::Dynamic
 * @param m specified matrix, vector, or row vector
 * @return maximum coefficient value in the vector, or infinity if is size zero
 */
template <typename T, int R, int C>
inline T max(const Eigen::Matrix<T, R, C>& m) {
  if (m.size() == 0) {
    return -std::numeric_limits<double>::infinity();
  }
  return m.maxCoeff();
}

}  // namespace math
}  // namespace stan
#endif
