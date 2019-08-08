#ifndef STAN_MATH_PRIM_MAT_FUN_INVERSE_HPP
#define STAN_MATH_PRIM_MAT_FUN_INVERSE_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/err/check_square.hpp>

namespace stan {
namespace math {

/**
 * Returns the inverse of the specified matrix.
 * @param m Specified matrix.
 * @return Inverse of the matrix.
 */
template <typename T, int R, int C>
inline Eigen::Matrix<T, R, C> inverse(const Eigen::Matrix<T, R, C>& m) {
  check_square("inverse", "m", m);
  return m.inverse();
}

template <typename Derived>
using scalar_t = typename Derived::Scalar;

template <typename Derived>
using eigen_return_t = typename Derived::PlainObject;

template <typename Derived>
using eigen_sparse_matrix = typename Eigen::SparseMatrix<scalar_t<Derived>>;

template <typename Derived>
using eigen_sparse_matrix_input = typename Eigen::SparseMatrixBase<Derived>;

template <typename Derived>
inline eigen_return_t<Derived> inverse(const eigen_sparse_matrix_input<Derived>& m) {
  Eigen::SimplicialLLT<eigen_sparse_matrix<Derived>> m_solver;
  m_solver.compute(m);
  eigen_sparse_matrix<Derived> identity(m.rows(), m.cols());
  identity.setIdentity();
  return m_solver.solve(identity);
}

}  // namespace math
}  // namespace stan
#endif
