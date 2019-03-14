#ifndef STAN_MATH_PRIM_FUN_DIVIDE_HPP
#define STAN_MATH_PRIM_FUN_DIVIDE_HPP

#include <stan/math/prim/err/domain_error.hpp>
#include <stan/math/prim/meta/likely.hpp>
#include <stan/math/prim/meta/return_type.hpp>
#include <cstddef>
#include <cstdlib>
#include <stan/math/prim/fun/Eigen.hpp>
#include <type_traits>











namespace stan {
namespace math {

/**
 * Return the division of the first scalar by
 * the second scalar.
 * @param[in] x Specified vector.
 * @param[in] y Specified scalar.
 * @return Vector divided by the scalar.
 */
template <typename T1, typename T2>
inline typename stan::return_type<T1, T2>::type divide(const T1& x,
                                                       const T2& y) {
  return x / y;
}

inline int divide(int x, int y) {
  if (unlikely(y == 0))
    domain_error("divide", "denominator is", y, "");
  return x / y;
}













/**
 * Return specified matrix divided by specified scalar.
 * @tparam R Row type for matrix.
 * @tparam C Column type for matrix.
 * @param m Matrix.
 * @param c Scalar.
 * @return Matrix divided by scalar.
 */
template <int R, int C, typename T>
inline typename std::enable_if<std::is_arithmetic<T>::value,
                               Eigen::Matrix<double, R, C> >::type
divide(const Eigen::Matrix<double, R, C>& m, T c) {
  return m / c;
}

}  // namespace math
}  // namespace stan
#endif
