#ifndef STAN_MATH_PRIM_FUN_HEAD_HPP
#define STAN_MATH_PRIM_FUN_HEAD_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Return the specified number of elements as a vector or row vector (same as
 * input) from the front of the specified vector or row vector.
 *
 * @tparam T type of the vector
 * @param v Vector input.
 * @param n Size of return.
 * @return The first n elements of v.
 * @throw std::out_of_range if n is out of range.
 */
template <typename T, require_vector_t<T>* = nullptr>
inline auto head(const T& v, size_t n) {
  if (n != 0) {
    check_vector_index("head", "n", v, n);
  }
  return v.head(n);
}

/**
 * Return the specified number of elements as a standard vector
 * from the front of the specified standard vector.
 *
 * @tparam T type of elements in the vector
 * @param sv Standard vector.
 * @param n Size of return.
 * @return The first n elements of sv.
 * @throw std::out_of_range if n is out of range.
 */
template <typename T>
std::vector<T> head(const std::vector<T>& sv, size_t n) {
  if (n != 0) {
    check_std_vector_index("head", "n", sv, n);
  }

  std::vector<T> s(sv.begin(), sv.begin() + n);
  return s;
}

}  // namespace math
}  // namespace stan

#endif
