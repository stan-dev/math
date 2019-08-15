#ifndef STAN_MATH_PRIM_MAT_FUN_CROSSPROD_HPP
#define STAN_MATH_PRIM_MAT_FUN_CROSSPROD_HPP

#include <stan/math/prim/mat/fun/typedefs.hpp>
#include <stan/math/prim/mat/fun/tcrossprod.hpp>

namespace stan {
namespace math {

/**
 * Returns the result of pre-multiplying a matrix by its
 * own transpose.
 * @param M Matrix to multiply.
 * @return Transpose of M times M
 */
template <typename T, enable_if_eigen<T>* = nullptr>
inline auto crossprod(const T& M) {
  return tcrossprod(M.transpose());
}

}  // namespace math
}  // namespace stan
#endif
