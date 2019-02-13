#ifndef STAN_MATH_PRIM_MAT_ERR_IS_CHOLESKY_FACTOR_CORR_HPP
#define STAN_MATH_PRIM_MAT_ERR_IS_CHOLESKY_FACTOR_CORR_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/scal/err/is_pos.hpp>
#include <stan/math/prim/mat/err/is_lower_triangular.hpp>
#include <stan/math/prim/mat/err/is_square.hpp>
#include <stan/math/prim/mat/err/constraint_tolerance.hpp>
#include <stan/math/prim/mat/err/is_unit_vector.hpp>

namespace stan {
namespace math {

/**
 * Check if the specified matrix is a
 * Cholesky factor of a correlation matrix.
 *
 * A Cholesky factor is a lower triangular matrix
 * whose diagonal elements are all positive. ThiS
 * definition does not require a square matrix just
 * that M &gt;= N, for M rows and N columns.
 *
 * Tolerance is specified by <code>math::CONSTRAINT_TOLERANCE</code>
 * as 1E-8.
 *
 * @tparam T_y Type of elements of Cholesky factor
 *
 * @param y Matrix to test
 *
 * @return <code>true</code>
 */
template <typename T_y>
inline bool is_cholesky_factor_corr(
    const Eigen::Matrix<T_y, Eigen::Dynamic, Eigen::Dynamic>& y) {
  using Eigen::Dynamic;
}

}  // namespace math
}  // namespace stan
#endif
