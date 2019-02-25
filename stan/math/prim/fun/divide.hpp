#ifndef STAN_MATH_PRIM_FUN_DIVIDE_HPP
#define STAN_MATH_PRIM_FUN_DIVIDE_HPP

#include <stanh/prim/err/domain_error.hpp>
#include <stanh/prim/meta/likely.hpp>
#include <stanh/prim/meta/return_type.hpp>
#include <cstddef>
#include <cstdlib>

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

}  // namespace math
}  // namespace stan
#endif
#ifndef STAN_MATH_PRIM_FUN_DIVIDE_HPP
#define STAN_MATH_PRIM_FUN_DIVIDE_HPP

#include <boost/type_traits/is_arithmetic.hpp>
#include <boost/utility/enable_if.hpp>
#include <stanh/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Return specified matrix divided by specified scalar.
 * @tparam R Row type for matrix.
 * @tparam C Column type for matrix.
 * @param m Matrix.
 * @param c Scalar.
 * @return Matrix divided by scalar.
 */
template <int R, int C, typename T>
inline typename boost::enable_if_c<boost::is_arithmetic<T>::value,
                                   Eigen::Matrix<double, R, C> >::type
divide(const Eigen::Matrix<double, R, C>& m, T c) {
  return m / c;
}

}  // namespace math
}  // namespace stan
#endif
#ifndef STAN_MATH_PRIM_FUN_DIVIDE_HPP
#define STAN_MATH_PRIM_FUN_DIVIDE_HPP

#include <stanh/prim/err/domain_error.hpp>
#include <stanh/prim/meta/likely.hpp>
#include <stanh/prim/meta/return_type.hpp>
#include <cstddef>
#include <cstdlib>

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

}  // namespace math
}  // namespace stan
#endif
