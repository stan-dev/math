#ifndef STAN_MATH_PRIM_FUN_IDENTITY_MATRIX_HPP
#define STAN_MATH_PRIM_FUN_IDENTITY_MATRIX_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Return a square identity matrix
 *
 * @param K size of the matrix
 * @return An identity matrix of size K.
 * @throw std::domain_error if K is negative.
 */
template <typename T = Eigen::MatrixXd, require_eigen_t<T>* = nullptr>
inline auto identity_matrix(int K) {
  check_nonnegative("identity_matrix", "size", K);
  return T::Identity(K, K);
}

}  // namespace math
}  // namespace stan

#endif
