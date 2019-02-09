#ifndef STAN_MATH_PRIM_MAT_FUN_CHOLESKY_DECOMPOSE_HPP
#define STAN_MATH_PRIM_MAT_FUN_CHOLESKY_DECOMPOSE_HPP

#include <stan/math/gpu/opencl_context.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/err/check_pos_definite.hpp>
#include <stan/math/prim/mat/err/check_square.hpp>
#include <stan/math/prim/mat/err/check_symmetric.hpp>
#ifdef STAN_OPENCL
#include <stan/math/gpu/cholesky_decompose.hpp>
#endif

#include <cmath>

namespace stan {
namespace math {

/**
 * Return the lower-triangular Cholesky factor (i.e., matrix
 * square root) of the specified square, symmetric matrix.  The return
 * value \f$L\f$ will be a lower-traingular matrix such that the
 * original matrix \f$A\f$ is given by
 * <p>\f$A = L \times L^T\f$.
 * @param m Symmetrix matrix.
 * @return Square root of matrix.
 * @throw std::domain_error if m is not a symmetric matrix or
 *   if m is not positive definite (if m has more than 0 elements)
 */
template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> cholesky_decompose(
    const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& m) {
  check_square("cholesky_decompose", "m", m);
  check_symmetric("cholesky_decompose", "m", m);

#ifdef STAN_OPENCL
  if (m.rows() >= opencl_context.tuning_opts().cholesky_size_worth_transfer) {
    std::cout  << "runs primitve cholesky without OpenCL" << std::endl;
    Eigen::LLT<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > llt(m.rows());
    matrix_gpu m_gpu(m);
    llt.compute(m);
    check_pos_definite("cholesky_decompose", "m", llt);
    
    
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> m_chol(m.rows(), m.cols());
    cholesky_decompose(m_gpu);
    copy(m_chol, m_gpu);  // NOLINT    
    for(int i=0;i<m.rows();i++){
      for(int j=0;j<=i;j++){
      std::cout << llt.matrixL()(i,j) << "\t" << m_chol(i,j) << std::endl;
      } 
    }
    return m_chol;
  } else {
    std::cout  << "runs primitve cholesky without OpenCL" << std::endl;
    Eigen::LLT<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > llt(m.rows());
    llt.compute(m);
    check_pos_definite("cholesky_decompose", "m", llt);
    return llt.matrixL();
  }
#else
  Eigen::LLT<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > llt(m.rows());
  llt.compute(m);
  check_pos_definite("cholesky_decompose", "m", llt);
  return llt.matrixL();
#endif
}
}  // namespace math

}  // namespace stan
#endif
