#ifndef STAN_MATH_PRIM_MAT_ERR_IS_MAT_VECTOR_HPP
#define STAN_MATH_PRIM_MAT_ERR_IS_MAT_VECTOR_HPP

#include <stan/math/prim/scal/meta/scalar_type.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Return <code>true</code> if the runtime size of the matrix **y**
 * is a row or column vector.
 * @tparam T Scalar type of the matrix, requires class method
 *   <code>.rows()</code> and <code>.cols()</code>
 * @tparam R Compile time rows of the matrix
 * @tparam C Compile time columns of the matrix
 * @param y Matrix to test
 * @return <code>true</code> if **y** is a row or column vector
 **/
template <typename T, int R, int C>
inline bool is_mat_vector(const Eigen::Matrix<T, R, C>& y) {
  return (R == 1) || (C == 1) || (y.rows() == 1 || y.cols() == 1);
}

}  // namespace math
}  // namespace stan
#endif
