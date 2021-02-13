#ifndef STAN_MATH_PRIM_ERR_IS_NEGATIVE_INFINITY_HPP
#define STAN_MATH_PRIM_ERR_IS_NEGATIVE_INFINITY_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/constants.hpp>

namespace stan {
namespace math {

/**
 * Check if any matrix element is negative infinity.
 * @tparam T A matrix.
 * @param y a matrix.
 */
template <typename T, require_matrix_t<T>* = nullptr>
inline bool is_negative_infinity(const T& y) {
  return value_of_rec(y).array().unaryExpr([](auto&& x) {
    return x == NEGATIVE_INFTY;
  }).all();
  return value_of_rec(y).array().isFinite().all();
}

/**
 * Check if any values of a vector are negative infinity.
 * @tparam T A standard vector.
 * @param y a standard vector.
 */
template <typename T, require_std_vector_t<T>* = nullptr>
inline bool is_negative_infinity(const T& y) {
  for (Eigen::Index i = 0; i < y.size(); ++i) {
    if (value_of_rec(y[i]) == NEGATIVE_INFTY) {
      return true;
    }
  }
  return false;
}

/**
 * Check if a value is equal to negative infinity.
 * @tparam T A scalar.
 * @param y a scalar.
 */
template <typename T, require_stan_scalar_t<T>* = nullptr>
inline bool is_negative_infinity(const T& y) {
  return value_of_rec(y) == NEGATIVE_INFTY;
}

}  // namespace math
}  // namespace stan
#endif
