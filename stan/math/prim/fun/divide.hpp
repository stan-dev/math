#ifndef STAN_MATH_PRIM_FUN_DIVIDE_HPP
#define STAN_MATH_PRIM_FUN_DIVIDE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <cstddef>
#include <cstdlib>
#include <type_traits>

namespace stan {
namespace math {

/**
 * Return the division of the first scalar by
 * the second scalar.
 * @param[in] x Specified scalar.
 * @param[in] y Specified scalar.
 * @return Scalar divided by the scalar.
 */
template <typename T1, typename T2,
          typename = require_all_stan_scalar_t<T1, T2>>
inline return_type_t<T1, T2> divide(const T1& x, const T2& y) {
  return x / y;
}

inline int divide(int x, int y) {
  if (unlikely(y == 0)) {
    throw_domain_error("divide", "denominator is", y, "");
  }
  return x / y;
}

/**
 * Return matrix divided by scalar.
 *
 * @tparam T1 type of the matrix or expression
 * @tparam T2 type of the scalar
 * @param[in] m specified matrix or expression
 * @param[in] c specified scalar
 * @return matrix divided by the scalar
 */
template <typename T1, typename T2, typename = require_eigen_t<T1>,
          typename = require_stan_scalar_t<T2>,
          typename = require_all_not_var_t<scalar_type_t<T1>, T2>>
inline auto divide(const T1& m, T2 c) {
  return (m / c).eval();
}

}  // namespace math
}  // namespace stan

#endif
