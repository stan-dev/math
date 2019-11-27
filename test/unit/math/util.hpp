#ifndef TEST_UNIT_MATH_UTIL_HPP
#define TEST_UNIT_MATH_UTIL_HPP

#include <stan/math.hpp>
#include <string>
#include <vector>

namespace stan {
namespace test {

/**
 * Return the Eigen vector with the same size and elements as the
 * specified standard vector.  Elements are copied from the specifed
 * input vector.
 *
 * @tparam T type of scalars in containers
 * @param x standard vector to copy
 * @return Eigen vector corresponding to specified standard vector
 */
template <typename T>
Eigen::Matrix<T, -1, 1> to_eigen_vector(const std::vector<T>& x) {
  return Eigen::Map<const Eigen::Matrix<T, -1, 1>>(x.data(), x.size());
}

/**
 * Return the standard vector with the same size and elements as the
 * specified Eigen matrix, vector, or row vector.
 *
 * @tparam T type of scalars in containers
 * @tparam R row specification of input matrix
 * @tparam C column specification of input matrix
 * @param x Eigen matrix, vector, or row vector to copy
 * @return standard vector corresponding to input
 */
template <typename T, int R, int C>
std::vector<T> to_std_vector(const Eigen::Matrix<T, R, C>& x) {
  std::vector<T> y;
  y.reserve(x.size());
  for (int i = 0; i < x.size(); ++i)
    y.push_back(x(i));
  return y;
}

}  // namespace test
}  // namespace stan
#endif
