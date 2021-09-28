#ifndef STAN_MATH_PRIM_FUN_FMA_HPP
#define STAN_MATH_PRIM_FUN_FMA_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the product of the first two arguments plus the third
 * argument.
 *
 * <p><i>Warning:</i> This does not delegate to the high-precision
 * platform-specific <code>fma()</code> implementation.
 *
 * @param x First argument.
 * @param y Second argument.
 * @param z Third argument.
 * @return The product of the first two arguments plus the third
 * argument.
 */
template <typename T1, typename T2, typename T3,
          require_all_arithmetic_t<T1, T2, T3>* = nullptr>
inline double fma(T1 x, T2 y, T3 z) {
  using std::fma;
  return fma(x, y, z);
}

template <typename T1, typename T2, typename T3,
          require_any_matrix_t<T1, T2, T3>* = nullptr,
          require_not_var_t<return_type_t<T1, T2, T3>>* = nullptr>
inline auto fma(T1&& x, T2&& y, T3&& z) {
  if (is_matrix<T1>::value && is_matrix<T2>::value) {
    check_matching_dims("fma", "x", x, "y", y);
  }
  if (is_matrix<T1>::value && is_matrix<T3>::value) {
    check_matching_dims("fma", "x", x, "z", z);
  } else if (is_matrix<T2>::value && is_matrix<T3>::value) {
    check_matching_dims("fma", "y", y, "z", z);
  }
  return make_holder(
      [](auto&& x, auto&& y, auto&& z) {
        return ((as_array_or_scalar(x) * as_array_or_scalar(y))
                + as_array_or_scalar(z))
            .matrix();
      },
      std::forward<T1>(x), std::forward<T2>(y), std::forward<T3>(z));
}

}  // namespace math
}  // namespace stan
#endif
