#ifndef STAN_MATH_PRIM_FUN_SYMMETRIZE_HPP
#define STAN_MATH_PRIM_FUN_SYMMETRIZE_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Return a symmetric matrix using elements from the lower triangular part of
 * the input matrix.
 *
 * @tparam T type of elements in the matrix
 * @param m Matrix.
 * @throw std:invalid_argument if the matrix is not square.
 */
template <typename T, require_eigen_t<T>* = nullptr>
inline Eigen::Matrix<value_type_t<T>, Eigen::Dynamic, Eigen::Dynamic>
symmetrize(const T& m) {
  check_square("symmetrize", "m", m);
  return m.template selfadjointView<Eigen::Lower>();
}

}  // namespace math
}  // namespace stan

#endif
