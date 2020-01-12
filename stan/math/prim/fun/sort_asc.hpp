#ifndef STAN_MATH_PRIM_FUN_SORT_ASC_HPP
#define STAN_MATH_PRIM_FUN_SORT_ASC_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <algorithm>
#include <vector>

namespace stan {
namespace math {

/**
 * Return the specified standard vector in ascending order.
 *
 * @tparam T Type of elements contained in vector.
 * @param xs Vector to order.
 * @return Vector in ascending order.
 * @throw std::domain_error If any of the values are NaN.
 */
template <typename T>
inline std::vector<T> sort_asc(std::vector<T> xs) {
  check_not_nan("sort_asc", "container argument", xs);
  std::sort(xs.begin(), xs.end());
  return xs;
}

/**
 * Return the specified vector in ascending order.
 *
 * @tparam T type of elements in the vector
 * @tparam R number of rows, can be Eigen::Dynamic
 * @tparam C number of columns, can be Eigen::Dynamic
 *
 * @param xs Vector to order.
 * @return Vector in ascending order.
 * @throw std::domain_error If any of the values are NaN.
 */
template <typename T, int R, int C>
inline Eigen::Matrix<T, R, C> sort_asc(Eigen::Matrix<T, R, C> xs) {
  check_not_nan("sort_asc", "container argument", xs);
  std::sort(xs.data(), xs.data() + xs.size());
  return xs;
}

}  // namespace math
}  // namespace stan

#endif
