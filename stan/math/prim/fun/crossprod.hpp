#ifndef STAN_MATH_PRIM_FUN_CROSSPROD_HPP
#define STAN_MATH_PRIM_FUN_CROSSPROD_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <stan/math/prim/fun/tcrossprod.hpp>

namespace stan {
namespace math {

/**
 * Returns the result of pre-multiplying a matrix by its
 * own transpose.
 *
 * @tparam T type of the matrix (must be derived from \c Eigen::MatrixBase)
 * @param M Matrix to multiply.
 * @return Transpose of M times M
 */
template <typename T, require_eigen_t<T>* = nullptr>
inline auto crossprod(const T& M) {
  return tcrossprod(M.transpose());
}

}  // namespace math
}  // namespace stan

#endif
