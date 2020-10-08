#ifndef STAN_MATH_PRIM_ERR_IS_MAT_FINITE_HPP
#define STAN_MATH_PRIM_ERR_IS_MAT_FINITE_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {

/**
 * Return <code>true</code> is the specified matrix is finite.
 * @tparam T Scalar type of the matrix, requires class method
 * <code>.allFinite()</code>
 * @tparam EigMat A type derived from `EigenBase`
 * @param y Matrix to test
 * @return <code>true</code> if the matrix is finite
 **/
template <typename EigMat, require_eigen_t<EigMat>* = nullptr>
inline bool is_mat_finite(const EigMat& y) {
  return y.array().isFinite().all();
}

}  // namespace math
}  // namespace stan
#endif
