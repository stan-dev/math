#ifndef STAN_MATH_REV_MAT_FUN_DIVIDE_HPP
#define STAN_MATH_REV_MAT_FUN_DIVIDE_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/mat/fun/to_var.hpp>
#include <stan/math/prim/meta.hpp>
namespace stan {
namespace math {

/**
 * Return the division of the specified column vector by
 * the specified scalar.
 * @param[in] v Specified vector.
 * @param[in] c Specified scalar.
 * @return Vector divided by the scalar.
 */
template <typename T1, typename T2, int R, int C,
          typename = require_any_var_t<T1, T2>>
inline Eigen::Matrix<var, R, C> divide(const Eigen::Matrix<T1, R, C>& v,
                                       const T2& c) {
  return to_var(v) / to_var(c);
}

}  // namespace math
}  // namespace stan
#endif
