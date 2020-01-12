#ifndef STAN_MATH_PRIM_FUN_TAIL_HPP
#define STAN_MATH_PRIM_FUN_TAIL_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Return the specified number of elements as a vector
 * from the back of the specified vector.
 *
 * @tparam T type of elements in the vector
 * @param v Vector input.
 * @param n Size of return.
 * @return The last n elements of v.
 * @throw std::out_of_range if n is out of range.
 */
template <typename T>
inline Eigen::Matrix<T, Eigen::Dynamic, 1> tail(
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& v, size_t n) {
  if (n != 0) {
    check_row_index("tail", "n", v, n);
  }
  return v.tail(n);
}

/**
 * Return the specified number of elements as a row vector
 * from the back of the specified row vector.
 *
 * @tparam T type of elements in the vector
 * @param rv Row vector.
 * @param n Size of return row vector.
 * @return The last n elements of rv.
 * @throw std::out_of_range if n is out of range.
 */
template <typename T>
inline Eigen::Matrix<T, 1, Eigen::Dynamic> tail(
    const Eigen::Matrix<T, 1, Eigen::Dynamic>& rv, size_t n) {
  if (n != 0) {
    check_column_index("tail", "n", rv, n);
  }
  return rv.tail(n);
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
  using idx_t = typename index_type<std::vector<T> >::type;
  if (n != 0) {
    check_std_vector_index("tail", "n", sv, n);
  }
  std::vector<T> s;
  for (idx_t i = sv.size() - n; i < sv.size(); ++i) {
    s.push_back(sv[i]);
  }
  return s;
}

}  // namespace math
}  // namespace stan

#endif
