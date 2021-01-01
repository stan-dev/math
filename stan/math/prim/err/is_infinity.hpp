#ifndef STAN_MATH_PRIM_ERR_IS_INFINITY_HPP
#define STAN_MATH_PRIM_ERR_IS_INFINITY_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/get.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return <code>true</code> if <code>y</code> is finite.
 * This function is vectorized and will check each element of
 * <code>y</code>.
 * @tparam T_y Type of y
 * @param y Variable to check
 * @throw <code>true</code> if y is not infinity, -infinity, or NaN
 */
template <typename T_y, require_matrix_t<T_y>* = nullptr>
inline bool is_infinity(const T_y& y) {
  return (y.array() == INFTY).any();
}

/**
 * Return <code>true</code> if <code>y</code> is finite.
 * This function is vectorized and will check each element of
 * <code>y</code>.
 * @tparam T_y Type of y
 * @param y Variable to check
 * @throw <code>true</code> if y is not infinity, -infinity, or NaN
 */
template <typename T_y, require_std_vector_t<T_y>* = nullptr>
inline bool is_infinity(const T_y& y) {
  return (as_array_or_scalar(y) == INFTY).any();
}

/**
 * Return <code>true</code> if <code>y</code> is finite.
 * This function is vectorized and will check each element of
 * <code>y</code>.
 * @tparam T_y Type of y
 * @param y Variable to check
 * @throw <code>true</code> if y is not infinity, -infinity, or NaN
 */
template <typename T_y, require_stan_scalar_t<T_y>* = nullptr>
inline bool is_infinity(const T_y& y) {
  return value_of_rec(y) == INFTY;
}

}  // namespace math
}  // namespace stan
#endif
