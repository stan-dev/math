#ifndef STAN_MATH_PRIM_SCAL_ERR_IS_CONSISTENT_SIZES_HPP
#define STAN_MATH_PRIM_SCAL_ERR_IS_CONSISTENT_SIZES_HPP

#include <stan/math/prim/scal/err/is_consistent_size.hpp>
#include <algorithm>

namespace stan {
namespace math {

/**
 * Check if the dimension of <code>x1</code> is consistent with
 * <code>x2</code>.
 * Consistent size is defined as having the same size if vector-like or
 * being a scalar.
 * @tparam T1 Type of <code>x1</code>
 * @tparam T2 Type of <code>x2</code>
 * @param x1 Variable to check for consistent size
 * @param name2 Variable name (for error messages)
 * @param x2 Variable to check for consistent size
 * @return <code>true</code> if sizes are consistent
 */
template <typename T1, typename T2>
inline bool is_consistent_sizes(const T1& x1, const char* name2, const T2& x2) {
  size_t max_size = std::max(is_vector<T1>::value * size_of(x1),
                             is_vector<T2>::value * size_of(x2));
  return (is_consistent_size(x1, max_size) && is_consistent_size(x2, max_size));
}

/**
 * Check if the dimension of <code>x1</code>, <code>x2</code>, and
 * <code>x3</code> are consistent.
 * Consistent size is defined as having the same size if vector-like or
 * being a scalar.
 * @tparam T1 Type of <code>x1</code>
 * @tparam T2 Type of <code>x2</code>
 * @tparam T3 Type of <code>x3</code>
 * @param x1 Variable to check for consistent size
 * @param x2 Variable to check for consistent size
 * @param x3 Variable to check for consistent size
 * @return <code>true</code> if sizes are consistent
 */
template <typename T1, typename T2, typename T3>
inline bool is_consistent_sizes(const T1& x1, const T2& x2, const T3& x3) {
  size_t max_size = std::max(is_vector<T1>::value * size_of(x1),
                             std::max(is_vector<T2>::value * size_of(x2),
                                      is_vector<T3>::value * size_of(x3)));
  return (is_consistent_size(x1, max_size) && is_consistent_size(x2, max_size)
          && is_consistent_size(x3, max_size));
}

/**
 * Check if the dimension of <code>x1</code>, <code>x2</code>, <code>x3</code>,
 * and <code>x4</code> are consistent.
 * Consistent size is defined as having the same size if vector-like or being a
 * scalar.
 * @tparam T1 Type of <code>x1</code>
 * @tparam T2 Type of <code>x2</code>
 * @tparam T3 Type of <code>x3</code>
 * @tparam T4 Type of <code>x4</code>
 * @param x1 Variable to check for consistent size
 * @param x2 Variable to check for consistent size
 * @param x3 Variable to check for consistent size
 * @param x4 Variable to check for consistent size
 * @return <code>true</code> if sizes are consistent
 */
template <typename T1, typename T2, typename T3, typename T4>
inline bool is_consistent_sizes(const T1& x1, const T2& x2, const T3& x3,
                                const T4& x4) {
  size_t max_size
      = std::max(is_vector<T1>::value * size_of(x1),
                 std::max(is_vector<T2>::value * size_of(x2),
                          std::max(is_vector<T3>::value * size_of(x3),
                                   is_vector<T4>::value * size_of(x4))));
  return (is_consistent_size(x1, max_size) && is_consistent_size(x2, max_size)
          && is_consistent_size(x3, max_size)
          && is_consistent_size(x4, max_size));
}

template <typename T1, typename T2, typename T3, typename T4, typename T5>
inline bool is_consistent_sizes(const char* function, const T1& x1,
                                const T2& x2, const T3& x3, const T4& x4,
                                const T5& x5) {
  size_t max_size = std::max(
      size_of(x1),
      std::max(size_of(x2),
               std::max(size_of(x3), std::max(size_of(x4), size_of(x5)))));
  return (is_consistent_size(x1, max_size) && is_consistent_size(x2, max_size)
          && is_consistent_size(x3, max_size)
          && is_consistent_size(x4, max_size)
          && is_consistent_size(x5, max_size));
}

}  // namespace math
}  // namespace stan
#endif
