#ifndef STAN_MATH_PRIM_MAT_ERR_IS_SIMPLEX_HPP
#define STAN_MATH_PRIM_MAT_ERR_IS_SIMPLEX_HPP

#include <stan/math/prim/arr/err/is_nonzero_size.hpp>
#include <stan/math/prim/mat/err/constraint_tolerance.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/meta/index_type.hpp>
#include <stan/math/prim/scal/meta/error_index.hpp>

namespace stan {
namespace math {

/**
 * Return <code>true</code> if the specified vector is simplex. To be a simplex,
 * all values must be greater than or equal to 0 and the values must sum to 1.
 * A valid simplex is one where the sum of hte elements is equal to 1. This
 * function tests that the sum is within the tolerance specified by
 * <code>CONSTRAINT_TOLERANCE</code>. This function only accepts Eigen vectors,
 * statically typed vectors, not general matrices with 1 column.
 * @tparam T_prob Scalar type of the vector
 * @param theta Vector to test
 * @return <code>true</code> if <code>theta</code> is not a 0-vector,
 *   if the vector is a simplex and if no element is <code>NaN</code>
 */
template <typename T_prob>
inline bool is_simplex(const Eigen::Matrix<T_prob, Eigen::Dynamic, 1>& theta) {
  return is_nonzero_size(theta)
    && fabs(1 - theta.sum()) < 1e-8
    && (theta.array() >= 0).all();
}

}  // namespace math
}  // namespace stan
#endif
