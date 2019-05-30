#ifndef STAN_MATH_PRIM_MAT_ERR_IS_SPSD_MATRIX_HPP
#define STAN_MATH_PRIM_MAT_ERR_IS_SPSD_MATRIX_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/scal/err/is_positive.hpp>
#include <stan/math/prim/mat/err/is_pos_semidefinite.hpp>
#include <stan/math/prim/mat/err/is_symmetric.hpp>
#include <stan/math/prim/mat/err/is_square.hpp>

namespace stan {
namespace math {
/**
 * Return <code>true</code> if the specified matrix is a square, symmetric, and
 * positive semi-definite.
 * @tparam T Scalar type of the matrix, requires class method
 * <code>.rows()</code>
 * @param y Matrix to test
 * @return <code>true</code> if the matrix is square or if the matrix is not 0x0
 * if the matrix is symmetric or if the matrix is positive semi-definite
 */
template <typename T_y>
inline bool is_spsd_matrix(
    const Eigen::Matrix<T_y, Eigen::Dynamic, Eigen::Dynamic>& y) {
  return is_positive(y.rows()) && is_symmetric(y) && is_pos_semidefinite(y);
}

}  // namespace math
}  // namespace stan
#endif
