#ifndef STAN_MATH_PRIM_ERR_CHECK_CONSISTENT_SIZES_HPP
#define STAN_MATH_PRIM_ERR_CHECK_CONSISTENT_SIZES_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/check_consistent_size.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <algorithm>

namespace stan {
namespace math {

/**
 * Check if the dimension of x1 is consistent with x2.
 * Consistent size is defined as having the same size if vector-like or
 * being a scalar.
 * @tparam T1 Type of x1
 * @tparam T2 Type of x2
 * @param function Function name (for error messages)
 * @param name1 Variable name (for error messages)
 * @param x1 Variable to check for consistent size
 * @param name2 Variable name (for error messages)
 * @param x2 Variable to check for consistent size
 * @throw <code>invalid_argument</code> if sizes are inconsistent
 */
template <typename T1, typename T2>
inline void check_consistent_sizes(const char* function, const char* name1,
                                   const T1& x1, const char* name2,
                                   const T2& x2) {
  size_t max_size = std::max(is_vector<T1>::value * stan::math::size(x1),
                             is_vector<T2>::value * stan::math::size(x2));
  check_consistent_size(function, name1, x1, max_size);
  check_consistent_size(function, name2, x2, max_size);
}

/**
 * Check if the dimension of x1, x2, and x3 are consistent.
 * Consistent size is defined as having the same size if vector-like or
 * being a scalar.
 * @tparam T1 Type of x1
 * @tparam T2 Type of x2
 * @tparam T3 Type of x3
 * @param function Function name (for error messages)
 * @param name1 Variable name (for error messages)
 * @param x1 Variable to check for consistent size
 * @param name2 Variable name (for error messages)
 * @param x2 Variable to check for consistent size
 * @param name3 Variable name (for error messages)
 * @param x3 Variable to check for consistent size
 * @throw <code>invalid_argument</code> if sizes are inconsistent
 */
template <typename T1, typename T2, typename T3>
inline void check_consistent_sizes(const char* function, const char* name1,
                                   const T1& x1, const char* name2,
                                   const T2& x2, const char* name3,
                                   const T3& x3) {
  size_t max_size
      = std::max(is_vector<T1>::value * stan::math::size(x1),
                 std::max(is_vector<T2>::value * stan::math::size(x2),
                          is_vector<T3>::value * stan::math::size(x3)));
  check_consistent_size(function, name1, x1, max_size);
  check_consistent_size(function, name2, x2, max_size);
  check_consistent_size(function, name3, x3, max_size);
}

/**
 * Check if the dimension of x1, x2, x3, and x4 are consistent.
 * Consistent size is defined as having the same size if
 * vector-like or being a scalar.
 * @tparam T1 Type of x1
 * @tparam T2 Type of x2
 * @tparam T3 Type of x3
 * @tparam T4 Type of x4
 * @param function Function name (for error messages)
 * @param name1 Variable name (for error messages)
 * @param x1 Variable to check for consistent size
 * @param name2 Variable name (for error messages)
 * @param x2 Variable to check for consistent size
 * @param name3 Variable name (for error messages)
 * @param x3 Variable to check for consistent size
 * @param name4 Variable name (for error messages)
 * @param x4 Variable to check for consistent size
 * @throw <code>invalid_argument</code> if sizes are inconsistent
 */
template <typename T1, typename T2, typename T3, typename T4>
inline void check_consistent_sizes(const char* function, const char* name1,
                                   const T1& x1, const char* name2,
                                   const T2& x2, const char* name3,
                                   const T3& x3, const char* name4,
                                   const T4& x4) {
  size_t max_size = std::max(
      is_vector<T1>::value * stan::math::size(x1),
      std::max(is_vector<T2>::value * stan::math::size(x2),
               std::max(is_vector<T3>::value * stan::math::size(x3),
                        is_vector<T4>::value * stan::math::size(x4))));
  check_consistent_size(function, name1, x1, max_size);
  check_consistent_size(function, name2, x2, max_size);
  check_consistent_size(function, name3, x3, max_size);
  check_consistent_size(function, name4, x4, max_size);
}
template <typename T1, typename T2, typename T3, typename T4, typename T5>
inline void check_consistent_sizes(const char* function, const char* name1,
                                   const T1& x1, const char* name2,
                                   const T2& x2, const char* name3,
                                   const T3& x3, const char* name4,
                                   const T4& x4, const char* name5,
                                   const T5& x5) {
  size_t max_size = std::max(
      stan::math::size(x1),
      std::max(stan::math::size(x2),
               std::max(stan::math::size(x3),
                        std::max(stan::math::size(x4), stan::math::size(x5)))));
  check_consistent_size(function, name1, x1, max_size);
  check_consistent_size(function, name2, x2, max_size);
  check_consistent_size(function, name3, x3, max_size);
  check_consistent_size(function, name4, x4, max_size);
  check_consistent_size(function, name5, x5, max_size);
}

}  // namespace math
}  // namespace stan
#endif
