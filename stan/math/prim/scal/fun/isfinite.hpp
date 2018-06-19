#ifndef STAN_MATH_PRIM_SCAL_FUN_ISFINITE_HPP
#define STAN_MATH_PRIM_SCAL_FUN_ISFINITE_HPP

#include <stan/math/prim/scal/fun/val.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Checks if the given number has finite value.
 *
 * Return <code>true</code> if the specified variable's
 * value is finite.
 *
 * Termination happens when a non-stan type is
 * returned from val, preventing ADL from finding
 * this declaration of isfinite again.
 *
 * Needed for Eigen's fullPivLu() method on
 * non-Hermitian  matrix types. Called on stan
 * types from Eigen::complex::isfinite_impl via ADL.
 *
 * @tparam T type of stan object
 * @param t Variable to test.
 * @return <code>true</code> if variable is finite.
 */
template <class T>
inline bool isfinite(T const& t) {
  return isfinite(val(t));
}

}  // namespace math
}  // namespace stan
#endif
