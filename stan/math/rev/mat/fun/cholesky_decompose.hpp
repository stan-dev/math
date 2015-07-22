#ifndef STAN_MATH_REV_MAT_FUN_CHOLESKY_DECOMPOSE_HPP
#define STAN_MATH_REV_MAT_FUN_CHOLESKY_DECOMPOSE_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/fun/cholesky_decompose.hpp>
#include <stan/math/rev/scal/fun/value_of_rec.hpp>
#include <stan/math/prim/mat/fun/value_of_rec.hpp>
#include <stan/math/prim/mat/err/check_pos_definite.hpp>
#include <stan/math/prim/mat/err/check_square.hpp>
#include <stan/math/prim/mat/err/check_symmetric.hpp>
#include <Eigen/KroneckerProduct>

namespace stan {
  namespace math {
      class cholesky_decompose_v_vari : public vari {
      public:
        int M_;  // A.rows() = A.cols() = B.rows()
        double* A_;
        double* L_;
        vari** _variRefA;
        vari** _variRefL;

        cholesky_decompose_v_vari(const Eigen::Matrix<var, -1, -1>& A,
																	const Eigen::Matrix<double, -1, -1>& L_A)
          : vari(0.0),
            M_(A.rows()),
            A_(reinterpret_cast<double*>
               (stan::math::ChainableStack::memalloc_
                .alloc(sizeof(double) * A.rows() * A.cols()))),
            L_(reinterpret_cast<double*>
               (stan::math::ChainableStack::memalloc_
                .alloc(sizeof(double) * A.rows() * (A.rows() + 1) / 2))),
            _variRefA(reinterpret_cast<vari**>
                      (stan::math::ChainableStack::memalloc_
                       .alloc(sizeof(vari*) * A.rows() * A.cols()))),
            _variRefL(reinterpret_cast<vari**>
                      (stan::math::ChainableStack::memalloc_
                       .alloc(sizeof(vari*) * A.rows() * (A.rows() + 1) / 2))) {
          using Eigen::Matrix;
          using Eigen::Map;

          size_t pos = 0;
          for (size_type j = 0; j < M_; j++) {
            for (size_type i = 0; i < M_; i++) {
              _variRefA[pos] = A(i, j).vi_;
              A_[pos++] = A(i, j).val();
            }
          }

          pos = 0;
          for (size_type j = 0; j < M_; j++) {
            for (size_type i = j; i < M_; i++) {
              L_[pos] = L_A(i, j);
              _variRefL[pos] = new vari(L_[pos], false);
              pos++;
            }
          }
        }

        virtual void chain() {
          using Eigen::Matrix;
          using Eigen::Map;
          Eigen::Matrix<double, R1, C1> adjA(M_, M_);
          Eigen::Matrix<double, R1, C2> adjL(M_, M_);

          size_t pos = 0;
          for (size_type j = 0; j < adjL.cols(); j++)
            for (size_type i = j; i < adjL.rows(); i++)
              adjL(i, j) = _variRefL[pos++]->adj_;

//          adjB = Map<Matrix<double, R1, C1> >(A_, M_, M_)
//            .transpose().colPivHouseholderQr().solve(adjC);
//          adjA.noalias() = -adjB
//            * Map<Matrix<double, R1, C2> >(C_, M_, N_).transpose();

          pos = 0;
          for (size_type j = 0; j < adjA.cols(); j++)
            for (size_type i = j; i < adjA.rows(); i++)
              _variRefA[pos++]->adj_ += adjA(i, j);
        }
      };

    Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic>
    cholesky_decompose(const Eigen::Matrix<var, -1, -1> &A) {
      Eigen::Matrix<var, -1, -1> L(A.rows(), A.cols());
			Eigen::Matrix<double, -1, -1> L_dbl(A.rows(), A.cols());
			
			res_dbl = cholesky_decompose(value_of_rec(A));
      // NOTE: this is not a memory leak, this vari is used in the
      // expression graph to evaluate the adjoint, but is not needed
      // for the returned matrix.  Memory will be cleaned up with the
      // arena allocator.

      cholesky_decompose_v_vari *baseVari
        = new cholesky_decompose_v_vari(A, L_A);

      size_t pos = 0;
      for (size_type j = 0; j < res.cols(); j++)
        for (size_type i = j; i < res.rows(); i++)
          res(i, j).vi_ = baseVari->_variRefL[pos++];

      return res;
    }
  }
}
#endif
