#ifndef STAN_MATH_PRIM_ERR_IS_SQUARE_HPP
#define STAN_MATH_PRIM_ERR_IS_SQUARE_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/is_size_match.hpp>

namespace stan {
namespace math {

/**
 * Return <code>true</code> if the matrix is square. This check allows 0x0
 * matrices.
 * @tparam EigMat A type derived from `EigenBase`
 * @param y Matrix to test
 * @return <code>true</code> if matrix is square
 */
template <typename EigMat, require_eigen_t<EigMat>* = nullptr>
inline bool is_square(const EigMat& y) {
  return is_size_match(y.rows(), y.cols());
}

}  // namespace math
}  // namespace stan
#endif
