#ifndef STAN_MATH_PRIM_ERR_CHECK_CONSISTENT_SIZES_MVT_HPP
#define STAN_MATH_PRIM_ERR_CHECK_CONSISTENT_SIZES_MVT_HPP

#include <stan/math/prim/err/check_consistent_size_mvt.hpp>
#include <stan/math/prim/fun/size_mvt.hpp>
#include <stan/math/prim/meta/require_generics.hpp>
#include <algorithm>

namespace stan {
namespace math {
namespace internal {

/** Base case of recursion, this function is a no-op. */
inline void check_consistent_with_size_mvt(const char*, size_t) { return; }

/**
 * Check that the provided inputs are either non-vectors, or vectors of the
 * specified size.
 * @tparam T1 type of first input
 * @tparam Ts type of other inputs
 * @param function function name (for error messages)
 * @param size the size that vector inputs must be
 * @param name1 name of variable corresponding to first input
 * @param x1 first input
 * @param names_and_xs more inputs
 * @throw `invalid_argument` if sizes are inconsistent
 */
template <typename T1, typename... Ts>
inline void check_consistent_with_size_mvt(const char* function, size_t size,
                                           const char* name1, const T1& x1,
                                           const Ts&... names_and_xs) {
  check_consistent_size_mvt(function, name1, x1, size);
  check_consistent_with_size_mvt(function, size, names_and_xs...);
}

/**
 * Compute the size of a vector, or produce 0 for a nonvector.
 * This overload handles the vector case.
 * @tparam T type of `x`
 * @param x vector
 * @return size of `x`
 */
template <typename T, typename = require_vector_t<T>>
inline size_t size_ignore_nonvector_mvt(const T& x) {
  return stan::math::size_mvt(x);
}

/**
 * This overload handles the nonvector case.
 * @tparam T type of non-vector argument
 * @return 0
 */
template <typename T, typename = require_not_vector_t<T>, typename = void>
inline size_t size_ignore_nonvector_mvt(const T&) {
  return 0;
}
}  // namespace internal

/**
 * Check if the dimensions of the inputs are consistent.
 * Consistent size is defined as having the same size if
 * vector-like or being a scalar.
 *
 * E.g.: check_consistent_sizes("some_function", "x1", x1, "x2", x2, etc.).
 *
 * @tparam T1 type of first input
 * @tparam Ts type of other inputs
 * @param function function name (for error messages)
 * @param name1 name of variable corresponding to first input
 * @param x1 first input
 * @param names_and_xs more inputs
 * @throw `invalid_argument` if sizes are inconsistent
 */
template <typename T1, typename... Ts>
inline void check_consistent_sizes_mvt(const char* function, const char* name1,
                                       const T1& x1,
                                       const Ts&... names_and_xs) {
  size_t max_size
      = std::max({internal::size_ignore_nonvector_mvt(x1),
                  internal::size_ignore_nonvector_mvt(names_and_xs)...});
  internal::check_consistent_with_size_mvt(function, max_size, name1, x1,
                                           names_and_xs...);
}

}  // namespace math
}  // namespace stan
#endif
