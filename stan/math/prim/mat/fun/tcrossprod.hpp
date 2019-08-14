#ifndef STAN_MATH_PRIM_MAT_FUN_TCROSSPROD_HPP
#define STAN_MATH_PRIM_MAT_FUN_TCROSSPROD_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/fun/typedefs.hpp>

namespace stan {
namespace math {

/**
 * Returns the result of post-multiplying a matrix by its
 * own transpose.
 * @param M Matrix to multiply.
 * @return M times its transpose.
 */
template <typename T, typename = enable_if_eigen<T>>
inline auto tcrossprod(const T& M) {
  return M * M.transpose();
}

}  // namespace math
}  // namespace stan
#endif
