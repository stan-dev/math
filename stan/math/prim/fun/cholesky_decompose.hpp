#ifndef STAN_MATH_PRIM_FUN_CHOLESKY_DECOMPOSE_HPP
#define STAN_MATH_PRIM_FUN_CHOLESKY_DECOMPOSE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#ifdef STAN_OPENCL
#include <stan/math/opencl/opencl.hpp>
#endif

#include <cmath>

namespace stan {
namespace math {

/**
 * Return the lower-triangular Cholesky factor (i.e., matrix
 * square root) of the specified square, symmetric matrix.  The return
 * value \f$L\f$ will be a lower-triangular matrix such that the
 * original matrix \f$A\f$ is given by
 * <p>\f$A = L \times L^T\f$.
 *
 * @tparam T type of elements in the matrix
 * @param m Symmetric matrix.
 * @return Square root of matrix.
 * @note Because OpenCL only works on doubles there are two
 * <code>cholesky_decompose</code> functions. One that works on doubles
 * and another that works on all other types (this one).
 * @throw std::domain_error if m is not a symmetric matrix or
 *   if m is not positive definite (if m has more than 0 elements)
 */
template <typename T, require_eigen_t<T>* = nullptr,
          require_not_same_vt<double, T>* = nullptr,
          require_not_eigen_vt<is_var,T>* = nullptr>
inline Eigen::Matrix<value_type_t<T>, T::RowsAtCompileTime,
                     T::ColsAtCompileTime>
cholesky_decompose(const T& m) {
  eval_return_type_t<T>& m_eval = m.eval();
  check_symmetric("cholesky_decompose", "m", m_eval);
  check_not_nan("cholesky_decompose", "m", m_eval);
  Eigen::LLT<Eigen::Matrix<value_type_t<T>, T::RowsAtCompileTime,
                           T::ColsAtCompileTime> >
      llt(m.rows());
  llt.compute(m_eval);
  check_pos_definite("cholesky_decompose", "m", llt);
  return llt.matrixL();
}

/**
 * Return the lower-triangular Cholesky factor (i.e., matrix
 * square root) of the specified square, symmetric matrix.  The return
 * value \f$L\f$ will be a lower-triangular matrix such that the
 * original matrix \f$A\f$ is given by
 * <p>\f$A = L \times L^T\f$.
 *
 * @param m Symmetric matrix.
 * @return Square root of matrix.
 * @note Because OpenCL only works on doubles there are two
 * <code>cholesky_decompose</code> functions. One that works on doubles
 * (this one) and another that works on all other types.
 * @throw std::domain_error if m is not a symmetric matrix or
 *   if m is not positive definite (if m has more than 0 elements)
 */
template <typename T, require_eigen_t<T>* = nullptr,
          require_same_vt<double, T>* = nullptr>
inline Eigen::Matrix<double, T::RowsAtCompileTime, T::ColsAtCompileTime>
cholesky_decompose(const T& m) {
  eval_return_type_t<T>& m_eval = m.eval();
  check_not_nan("cholesky_decompose", "m", m);
#ifdef STAN_OPENCL
  if (m.rows() >= opencl_context.tuning_opts().cholesky_size_worth_transfer) {
    matrix_cl<double> m_cl(m);
    return from_matrix_cl(cholesky_decompose(m_cl));
  } else {
    check_symmetric("cholesky_decompose", "m", m);
    Eigen::LLT<
        Eigen::Matrix<double, T::RowsAtCompileTime, T::ColsAtCompileTime> >
        llt(m.rows());
    llt.compute(m);
    check_pos_definite("cholesky_decompose", "m", llt);
    return llt.matrixL();
  }
#else
  check_symmetric("cholesky_decompose", "m", m);
  Eigen::LLT<Eigen::Matrix<double, T::RowsArCompileTime, T::ColsAtCompileTime> > llt(
      m.rows());
  llt.compute(m);
  check_pos_definite("cholesky_decompose", "m", llt);
  return llt.matrixL();
#endif
}

}  // namespace math
}  // namespace stan

#endif
