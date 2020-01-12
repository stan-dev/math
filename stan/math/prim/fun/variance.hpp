#ifndef STAN_MATH_PRIM_FUN_VARIANCE_HPP
#define STAN_MATH_PRIM_FUN_VARIANCE_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/mean.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Returns the sample variance (divide by length - 1) of the
 * coefficients in the specified standard vector.
 *
 * @tparam T type of elements in the vector
 * @param v specified vector
 * @return sample variance of vector
 * @throw <code>std::invalid_argument</code> if the vector has size zero
 */
template <typename T>
inline return_type_t<T> variance(const std::vector<T>& v) {
  check_nonzero_size("variance", "v", v);
  if (v.size() == 1) {
    return 0.0;
  }
  T v_mean(mean(v));
  T sum_sq_diff(0);
  for (size_t i = 0; i < v.size(); ++i) {
    T diff = v[i] - v_mean;
    sum_sq_diff += diff * diff;
  }
  return sum_sq_diff / (v.size() - 1);
}

/**
 * Returns the sample variance (divide by length - 1) of the
 * coefficients in the specified matrix
 *
 * @tparam T type of elements in the vector
 * @tparam R number of rows in the matrix, can be Eigen::Dynamic
 * @tparam C number of columns in the matrix, can be Eigen::Dynamic
 *
 * @param m matrix
 * @return sample variance of coefficients
 * @throw <code>std::invalid_argument</code> if the matrix has size zero
 */
template <typename T, int R, int C>
inline return_type_t<T> variance(const Eigen::Matrix<T, R, C>& m) {
  check_nonzero_size("variance", "m", m);

  if (m.size() == 1) {
    return 0.0;
  }
  return_type_t<T> mn(mean(m));
  return_type_t<T> sum_sq_diff(0);
  for (int i = 0; i < m.size(); ++i) {
    return_type_t<T> diff = m(i) - mn;
    sum_sq_diff += diff * diff;
  }
  return sum_sq_diff / (m.size() - 1);
}

}  // namespace math
}  // namespace stan

#endif
