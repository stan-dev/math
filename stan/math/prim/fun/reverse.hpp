#ifndef STAN_MATH_PRIM_FUN_REVERSE_HPP
#define STAN_MATH_PRIM_FUN_REVERSE_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <algorithm>
#include <vector>

namespace stan {
namespace math {

/**
 * Return a copy of the specified array in reversed order.
 *
 * @tparam T type of elements in the array
 * @param x array to reverse
 * @return Array in reversed order.
 */
template <typename T>
inline std::vector<T> reverse(const std::vector<T>& x) {
  std::vector<T> rev(x.size());
  std::reverse_copy(x.begin(), x.end(), rev.begin());
  return rev;
}

/**
 * Return a copy of the specified vector in reversed order.
 *
 * @tparam T type of elements in the vector
 * @param x vector to reverse
 * @return Vector in reversed order.
 */
template <typename T>
inline Eigen::Matrix<T, Eigen::Dynamic, 1> reverse(
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& x) {
  return x.reverse();
}

/**
 * Return a copy of the specified row vector in reversed order.
 *
 * @tparam T type of elements in the row vector
 * @param x row vector to reverse
 * @return Row vector in reversed order.
 */
template <typename T>
inline Eigen::Matrix<T, 1, Eigen::Dynamic> reverse(
    const Eigen::Matrix<T, 1, Eigen::Dynamic>& x) {
  return x.reverse();
}

}  // namespace math
}  // namespace stan
#endif
