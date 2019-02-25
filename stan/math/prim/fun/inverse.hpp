#ifndef STAN_MATH_PRIM_FUN_INVERSE_HPP
#define STAN_MATH_PRIM_FUN_INVERSE_HPP

#include <stanh/prim/fun/Eigen.hpp>
#include <stanh/prim/err/check_square.hpp>

namespace stan {
namespace math {

/**
 * Returns the inverse of the specified matrix.
 * @param m Specified matrix.
 * @return Inverse of the matrix.
 */
template <typename T, int R, int C>
inline Eigen::Matrix<T, R, C> inverse(const Eigen::Matrix<T, R, C>& m) {
  check_square("inverse", "m", m);
  return m.inverse();
}

}  // namespace math
}  // namespace stan
#endif
#ifndef STAN_MATH_PRIM_FUN_INVERSE_HPP
#define STAN_MATH_PRIM_FUN_INVERSE_HPP

#include <stanh/prim/fun/Eigen.hpp>
#include <stanh/prim/err/check_square.hpp>

namespace stan {
namespace math {

/**
 * Returns the inverse of the specified matrix.
 * @param m Specified matrix.
 * @return Inverse of the matrix.
 */
template <typename T, int R, int C>
inline Eigen::Matrix<T, R, C> inverse(const Eigen::Matrix<T, R, C>& m) {
  check_square("inverse", "m", m);
  return m.inverse();
}

}  // namespace math
}  // namespace stan
#endif
