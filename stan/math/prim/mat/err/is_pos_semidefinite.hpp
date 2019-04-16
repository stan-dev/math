#ifndef STAN_MATH_PRIM_MAT_ERR_IS_POS_SEMIDEFINITE_HPP
#define STAN_MATH_PRIM_MAT_ERR_IS_POS_SEMIDEFINITE_HPP

#include <stan/math/prim/mat/err/is_symmetric.hpp>
#include <stan/math/prim/mat/err/constraint_tolerance.hpp>
#include <stan/math/prim/scal/err/is_not_nan.hpp>
#include <stan/math/prim/scal/err/is_positive.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/meta/index_type.hpp>
#include <stan/math/prim/mat/fun/value_of_rec.hpp>

namespace stan {
namespace math {

/**
 * Return <code>true</code> if the specified matrix is positive definite.
 * @tparam T_y scalar type of the matrix, requires class method
 *<code>.rows()</code>
 * @param y Matrix to test
 * @return <code>true</code> if the matrix is square, or if the matrix does
 *   not have zero size, if the matrix is symmetric, or if it is positive
 *   semi-definite, or if no element of the matrix is <code>NaN</code>.
 **/
template <typename T_y>
inline bool is_pos_semidefinite(
    const Eigen::Matrix<T_y, Eigen::Dynamic, Eigen::Dynamic>& y) {
  if (!is_symmetric(y))
    return false;

  if (!is_positive(y.rows()))
    return false;

  if (y.rows() == 1 && !(y(0, 0) >= 0.0))
    return false;

  using Eigen::Dynamic;
  using Eigen::LDLT;
  using Eigen::Matrix;
  LDLT<Matrix<double, Dynamic, Dynamic> > cholesky = value_of_rec(y).ldlt();
  if (cholesky.info() != Eigen::Success
      || (cholesky.vectorD().array() < 0.0).any())
    return false;
  return is_not_nan(y);
}

}  // namespace math
}  // namespace stan
#endif
