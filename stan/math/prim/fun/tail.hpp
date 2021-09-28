#ifndef STAN_MATH_PRIM_FUN_TAIL_HPP
#define STAN_MATH_PRIM_FUN_TAIL_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Return the specified number of elements as a vector or row vector (same as
 * input) from the back of the specified vector or row vector.
 *
 * @tparam T type of the vector
 * @param v Vector input.
 * @param n Size of return.
 * @return The last n elements of v.
 * @throw std::out_of_range if n is out of range.
 */
template <typename T, require_vector_t<T>* = nullptr>
inline auto tail(const T& v, size_t n) {
  if (n != 0) {
    check_vector_index("tail", "n", v, n);
  }
  return v.tail(n);
}

/**
 * Return the specified number of elements as a standard vector
 * from the back of the specified standard vector.
 *
 * @tparam T type of elements in the vector
 * @param sv Standard vector.
 * @param n Size of return.
 * @return The last n elements of sv.
 * @throw std::out_of_range if n is out of range.
 */
template <typename T>
std::vector<T> tail(const std::vector<T>& sv, size_t n) {
  if (n != 0) {
    check_std_vector_index("tail", "n", sv, n);
  }
  std::vector<T> s(sv.end() - n, sv.end());
  return s;
}

}  // namespace math
}  // namespace stan

#endif
