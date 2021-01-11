#ifndef STAN_MATH_PRIM_ERR_IS_LDLT_FACTOR_HPP
#define STAN_MATH_PRIM_ERR_IS_LDLT_FACTOR_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/LDLT_factor.hpp>

namespace stan {
namespace math {

/**
 * Return <code>true</code> if the specified LDLT factor is invalid.
 * An <code>LDLT_factor</code> is invalid if it was constructed from
 * a matrix that is not positive definite.
 *
 * @tparam T Type matrix of LDLT
 * @param A The LDLT factor to check for validity
 * @return <code>true</code> if the LDLT factor is valid
 */
template <typename T>
inline bool is_ldlt_factor(LDLT_factor<T>& A) {
  return A.ldlt().info() == Eigen::Success && A.ldlt().isPositive()
         && (A.ldlt().vectorD().array() > 0).all();
}

}  // namespace math
}  // namespace stan
#endif
