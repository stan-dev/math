#ifndef STAN_MATH_PRIM_FUN_REVERSE_HPP
#define STAN_MATH_PRIM_FUN_REVERSE_HPP

#include <stan/math/prim/meta.hpp>
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
 * Return a copy of the specified vector or row vector
 * in reversed order.
 *
 * @tparam T type of container (vector or row vector)
 * @param x vector or row vector to reverse
 * @return Vector or row vector in reversed order.
 */
template <typename T, require_vector_t<T>* = nullptr,
  require_not_std_vector_t<T>* = nullptr>
inline auto reverse(T&& x) {
  return make_holder([](auto&& x_) {
    return x_.reverse();
  }, std::forward<T>(x));
}

}  // namespace math
}  // namespace stan
#endif
