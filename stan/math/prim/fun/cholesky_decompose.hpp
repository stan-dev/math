#ifndef STAN_MATH_PRIM_FUN_CHOLESKY_DECOMPOSE_HPP
#define STAN_MATH_PRIM_FUN_CHOLESKY_DECOMPOSE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

#include <cmath>

namespace stan {
namespace math {

/**
 * Return the lower-triangular Cholesky factor (i.e., matrix
 * square root) of the specified square, symmetric matrix.  The return
 * value \f$L\f$ will be a lower-triangular matrix such that the
 * original matrix \f$A\f$ is given by
 * <p>\f$A = L \times L^EigMat\f$.
 *
 * @tparam EigMat type of the matrix (must be derived from \c Eigen::MatrixBase)
 * @param m Symmetric matrix.
 * @return Square root of matrix.
 * @note Because OpenCL only works on doubles there are two
 * <code>cholesky_decompose</code> functions. One that works on doubles
 * and another that works on all other types (this one).
 * @throw std::domain_error if m is not a symmetric matrix or
 *   if m is not positive definite (if m has more than 0 elements)
 */
template <typename EigMat, require_eigen_t<EigMat>* = nullptr,
          require_not_eigen_vt<is_var, EigMat>* = nullptr>
inline Eigen::Matrix<value_type_t<EigMat>, EigMat::RowsAtCompileTime,
                     EigMat::ColsAtCompileTime>
cholesky_decompose(const EigMat& m) {
  const eval_return_type_t<EigMat>& m_eval = m.eval();
  check_symmetric("cholesky_decompose", "m", m_eval);
  check_not_nan("cholesky_decompose", "m", m_eval);
  Eigen::LLT<Eigen::Matrix<value_type_t<EigMat>, EigMat::RowsAtCompileTime,
                           EigMat::ColsAtCompileTime>>
      llt = m_eval.llt();
  check_pos_definite("cholesky_decompose", "m", llt);
  return llt.matrixL();
}

/**
 * Return the sparse lower-triangular Cholesky factor (i.e., matrix
 * square root) of the specified sparse square, symmetric matrix.  The return
 * value \f$L\f$ will be a sparse lower-triangular matrix such that the
 * original matrix \f$A\f$ is given by
 * <p>\f$A = P^TL  L^TP\f$,
 * where \f$P\f$ is a permutation matrix that has been computed to minimize
 * fill-in.
 *
 * @tparam SpEigMat type of the matrix (must be derived from \c
 * Eigen::SparseMatrixBase)
 * @param m Sparse symmetric matrix.
 * @return A Eigen::SimplicialLLT<SpEigMat> object containing
 * the sparse cholesky factor matrix and its associated metadata
 * (most importantly its permutation vector).
 * @note Unlike the dense matrix specialization, this version does
 * *not* return a matrix. The matrix can be accessed through the `matrixL()`
 * member function but remember this is the Cholesky factor of
 * @throw std::domain_error if m is not a symmetric matrix or
 *   if m is not positive definite (if m has more than 0 elements)
 */
template <typename SpEigMat, require_eigen_sparse_base_t<SpEigMat>* = nullptr,
          require_not_eigen_vt<is_var, SpEigMat>* = nullptr>
inline Eigen::SimplicialLLT<SpEigMat> cholesky_decompose(const EigMat& m) {
  const eval_return_type_t<SpEigMat>& m_eval = m.eval();
  check_symmetric("cholesky_decompose", "m", m_eval);
  check_not_nan("cholesky_decompose", "m", m_eval);
  Eigen::SimplicialLLT<SpEigMat> llt = m_eval.llt();
  check_pos_definite("cholesky_decompose", "m", llt);
  return llt;
}

}  // namespace math
}  // namespace stan

#endif
