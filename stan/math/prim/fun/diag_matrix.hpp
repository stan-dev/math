#ifndef STAN_MATH_PRIM_FUN_DIAG_MATRIX_HPP
#define STAN_MATH_PRIM_FUN_DIAG_MATRIX_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Return a square diagonal matrix with the specified vector of
 * coefficients as the diagonal values.
 *
 * @tparam EigVec type of the vector (must be derived from \c Eigen::MatrixBase
 * and have one compile time dimmension equal to 1)
 * @param[in] v Specified vector.
 * @return Diagonal matrix with vector as diagonal values.
 */
template <typename EigVec, require_eigen_vector_t<EigVec>* = nullptr>
inline Eigen::Matrix<value_type_t<EigVec>, Eigen::Dynamic, Eigen::Dynamic>
diag_matrix(const EigVec& v) {
  return v.asDiagonal();
}

}  // namespace math
}  // namespace stan

#endif
