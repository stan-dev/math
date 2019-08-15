#ifndef STAN_MATH_PRIM_MAT_FUN_INVERSE_HPP
#define STAN_MATH_PRIM_MAT_FUN_INVERSE_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/mat/err/check_square.hpp>

namespace stan {
namespace math {

/**
 * Returns the inverse of the specified matrix.
 * @param m Specified matrix.
 * @return Inverse of the matrix.
 */
template <typename T, enable_if_eigen<T>* = nullptr>
inline auto inverse(const T& m) {
  check_square("inverse", "m", m);
  return m.inverse();
}

}  // namespace math
}  // namespace stan
#endif
