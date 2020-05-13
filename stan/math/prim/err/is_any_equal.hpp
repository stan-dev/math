#ifndef STAN_MATH_PRIM_ERR_IS_ANY_EQUAL_HPP
#define STAN_MATH_PRIM_ERR_IS_ANY_EQUAL_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Return `true` if `x` is equal to `y`
 * @tparam T1 An Arithmetic type
 * @tparam T2 An Arithmetic type
 * @param x An Arithmetic value
 * @param y An Arithmetic value
 */
template <typename T1, typename T2, require_all_arithmetic_t<T1, T2>* = nullptr>
bool is_any_equal(T1 x, T2 y) {
  return x == y;
}

/**
 * Return `true` if any of `x`'s values are equal to `y`
 * @tparam T1 A type inheriting from `Eigen::EigenBase`
 * @tparam T2 An Arithmetic type
 * @param x An Eigen expression
 * @param y An Arithmetic value
 */
template <typename T1, typename T2, require_eigen_t<T1>* = nullptr,
          require_arithmetic_t<T2>* = nullptr>
bool is_any_equal(const T1& x, T2 y) {
  return (x.array() == y).any();
}

}  // namespace math
}  // namespace stan
#endif
