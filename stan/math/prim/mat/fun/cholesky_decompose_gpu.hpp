#ifndef STAN_MATH_PRIM_MAT_FUN_CHOLESKY_DECOMPOSE_GPU_HPP
#define STAN_MATH_PRIM_MAT_FUN_CHOLESKY_DECOMPOSE_GPU_HPP

#ifdef STAN_GPU
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/fun/ViennaCL.hpp>
#include <stan/math/prim/mat/err/check_pos_definite.hpp>
#include <stan/math/prim/mat/err/check_square.hpp>
#include <stan/math/prim/mat/err/check_symmetric.hpp>
#include <boost/type_traits/is_arithmetic.hpp>
#include <boost/utility/enable_if.hpp>
#include <stan/math/prim/scal/meta/contains_fvar.hpp>
#include <algorithm>

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
    typename boost::disable_if_c<stan::contains_fvar<T>::value || 
      stan::is_fvar<T>::value,
      Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >::type
    cholesky_decompose_gpu(const Eigen::Matrix
                       <T, Eigen::Dynamic, Eigen::Dynamic>& m) {
      check_square("cholesky_decompose", "m", m);
      check_symmetric("cholesky_decompose", "m", m);
      viennacl::matrix<T>  vcl_m(m.rows(), m.cols());
      viennacl::copy(m, vcl_m);
      viennacl::linalg::lu_factorize(vcl_m);
      Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> m_l(m.rows(), m.cols());
      viennacl::copy(vcl_m, m_l);
      // TODO(Steve/Sean): Where should this check go?
      // check_pos_definite("cholesky_decompose", "m", L_A);
      // STEVE LOOK HERE, NEED TO REMOVE DOUBLE 
      m_l = m_l.template triangularView<Eigen::Upper>().transpose();
      for (int i = 0; i < m_l.rows(); i++) m_l.col(i) /= std::sqrt(m_l(i, i));
      return m_l;
      //NOTE: (Steve/Sean) we need a check for positive definite in this call
    }

  }
}
#endif
#endif
