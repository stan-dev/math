#ifndef STAN_MATH_PRIM_ERR_CHECK_CONSISTENT_SIZES_MVT_HPP
#define STAN_MATH_PRIM_ERR_CHECK_CONSISTENT_SIZES_MVT_HPP

#include <stan/math/prim/err/check_consistent_size_mvt.hpp>
#include <stan/math/prim/fun/size_mvt.hpp>
#include <stan/math/prim/meta/require_generics.hpp>
#include <algorithm>

namespace stan {
namespace math {
namespace internal {
void check_consistent_sizes_mvt_impl(const char* function, size_t max_size) {
  return;
}

template <typename T1, typename... Ts>
void check_consistent_sizes_mvt_impl(const char* function, size_t max_size,
                                     const char* name1, const T1& x1,
                                     const Ts&... xs) {
  check_consistent_size_mvt(function, name1, x1, max_size);
  check_consistent_sizes_mvt_impl(function, max_size, xs...);
}

template <typename T, typename = require_vector_t<T>>
size_t size_mvt_ignore_nonvector(const T& x) {
  return stan::math::size_mvt(x);
}

template <typename T, typename = require_not_vector_t<T>, typename = void>
size_t size_mvt_ignore_nonvector(const T&) {
  return 0;
}
}  // namespace internal

/**
 * Check if the dimensions of the inputs are consistent.
 * Consistent size is defined as having the same size if
 * vector-like or being a scalar.
 *
 * E.g.: check_consistent_sizes_mvt("some_function", "x1", x1, "x2", x2, etc.).
 *
 * @tparam T1 type of first input
 * @tparam Ts type of other inputs
 * @param function Function name (for error messages)
 * @param name1 name of variable corresponding to first input
 * @param x1 first input
 * @param names_and_xs more inputs
 * @throw <code>invalid_argument</code> if sizes are inconsistent
 */
template <typename T1, typename... Ts>
inline void check_consistent_sizes_mvt(const char* function, const char* name1,
                                       const T1& x1,
                                       const Ts&... names_and_xs) {
  size_t max_size_mvt
      = std::max({internal::size_mvt_ignore_nonvector(x1),
                  internal::size_mvt_ignore_nonvector(names_and_xs)...});
  internal::check_consistent_sizes_mvt_impl(function, max_size_mvt, name1, x1,
                                            names_and_xs...);
}

}  // namespace math
}  // namespace stan
#endif
