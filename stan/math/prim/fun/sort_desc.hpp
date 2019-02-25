#ifndef STAN_MATH_PRIM_FUN_SORT_DESC_HPP
#define STAN_MATH_PRIM_FUN_SORT_DESC_HPP

#include <stanh/prim/fun/Eigen.hpp>
#include <stanh/prim/fun/value_of_rec.hpp>
#include <stanh/prim/err/check_not_nan.hpp>
#include <algorithm>
#include <functional>

namespace stan {
namespace math {

/**
 * Return the specified vector in descending order.
 *
 * @tparam T Type of elements contained in vector.
 * @param xs Vector to order.
 * @return Vector in descending order.
 * @throw std::domain_error If any of the values are NaN.
 */
template <typename T, int R, int C>
inline Eigen::Matrix<T, R, C> sort_desc(Eigen::Matrix<T, R, C> xs) {
  check_not_nan("sort_asc", "container argument", xs);
  std::sort(xs.data(), xs.data() + xs.size(), std::greater<T>());
  return xs;
}

}  // namespace math
}  // namespace stan
#endif
#ifndef STAN_MATH_PRIM_FUN_SORT_DESC_HPP
#define STAN_MATH_PRIM_FUN_SORT_DESC_HPP

#include <stanh/prim/err/check_not_nan.hpp>
#include <algorithm>
#include <functional>
#include <vector>

namespace stan {
namespace math {

/**
 * Return the specified standard vector in descending order.
 *
 * @tparam T Type of elements contained in vector.
 * @param xs Vector to order.
 * @return Vector in descending order.
 * @throw std::domain_error If any of the values are NaN.
 */
template <typename T>
inline std::vector<T> sort_desc(std::vector<T> xs) {
  check_not_nan("sort_asc", "container argument", xs);
  std::sort(xs.begin(), xs.end(), std::greater<T>());
  return xs;
}

}  // namespace math
}  // namespace stan
#endif
#ifndef STAN_MATH_PRIM_FUN_SORT_DESC_HPP
#define STAN_MATH_PRIM_FUN_SORT_DESC_HPP

#include <stanh/prim/fun/Eigen.hpp>
#include <stanh/prim/fun/value_of_rec.hpp>
#include <stanh/prim/err/check_not_nan.hpp>
#include <algorithm>
#include <functional>

namespace stan {
namespace math {

/**
 * Return the specified vector in descending order.
 *
 * @tparam T Type of elements contained in vector.
 * @param xs Vector to order.
 * @return Vector in descending order.
 * @throw std::domain_error If any of the values are NaN.
 */
template <typename T, int R, int C>
inline Eigen::Matrix<T, R, C> sort_desc(Eigen::Matrix<T, R, C> xs) {
  check_not_nan("sort_asc", "container argument", xs);
  std::sort(xs.data(), xs.data() + xs.size(), std::greater<T>());
  return xs;
}

}  // namespace math
}  // namespace stan
#endif
