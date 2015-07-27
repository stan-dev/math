#ifndef STAN_MATH_REV_MAT_FUN_CHOLESKY_DECOMPOSE_HPP
#define STAN_MATH_REV_MAT_FUN_CHOLESKY_DECOMPOSE_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/fun/cholesky_decompose.hpp>
#include <stan/math/rev/scal/fun/value_of_rec.hpp>
#include <stan/math/prim/mat/fun/value_of_rec.hpp>
#include <stan/math/prim/mat/err/check_pos_definite.hpp>
#include <stan/math/prim/mat/err/check_square.hpp>
#include <stan/math/prim/mat/err/check_symmetric.hpp>

namespace stan {
  namespace math {

    size_t indexer(int i, int j, int N) {
      return static_cast<size_t>(j * N + i);
    }

    class cholesky_decompose_v_vari : public vari {
    public:
      int M_;  // A.rows() = A.cols()
      double* L_;
      double* A_;
      vari** _variRefA;
      vari** _variRefL;
      vari** _dummy;

      cholesky_decompose_v_vari(const Eigen::Matrix<var, -1, -1>& A)
        : vari(0.0),
          M_(A.rows()),
          L_(reinterpret_cast<double*>
             (stan::math::ChainableStack::memalloc_
              .alloc(sizeof(double) * A.rows() * (A.rows() + 1) / 2))),
          A_(reinterpret_cast<double*>
             (stan::math::ChainableStack::memalloc_
              .alloc(sizeof(double) * A.rows() * A.rows()))),
          _variRefA(reinterpret_cast<vari**>
                    (stan::math::ChainableStack::memalloc_
                     .alloc(sizeof(vari*) * A.rows() * (A.cols() + 1) / 2))),
          _variRefL(reinterpret_cast<vari**>
                    (stan::math::ChainableStack::memalloc_
                     .alloc(sizeof(vari*) * A.rows() * (A.rows() + 1) / 2))),
          _dummy(reinterpret_cast<vari**>
                    (stan::math::ChainableStack::memalloc_
                     .alloc(sizeof(vari*))))  {
          using Eigen::Map;
          using Eigen::Matrix;

        size_t posA = 0;
        for (size_type j = 0; j < M_; ++j) 
          for (size_type k = 0; k < M_; ++k)  
            A_[posA++] = A.coeffRef(k, j).vi_->val_;

        Matrix<double, -1, -1> C = Map<Matrix<double, -1, -1> > (A_,M_,M_);
        Eigen::LLT<Eigen::MatrixXd> D = C.selfadjointView<Eigen::Lower>().llt();
        check_pos_definite("cholesky_decompose", "m", D);

        size_t pos = 0;
        for (size_type j = 0; j < M_; ++j) {
          for (size_type i = j; i < M_; ++i) {
            _variRefA[pos] = A.coeffRef(i, j).vi_;
            L_[pos] = D.matrixL()(i, j);
            _variRefL[pos] = new vari(L_[pos], false);
            pos++;
          }
        }
        _dummy[0] = new vari(0.0, false);
      }

      virtual void chain() {
        Eigen::Matrix<double, Dynamic, Dynamic> adjA(M_, M_);
        Eigen::Matrix<double, Dynamic, Dynamic> adjL(M_, M_);
        Eigen::Matrix<double, Dynamic, Dynamic> LA(M_, M_);
        size_t pos = -1;
        for (size_type j = 0; j < M_; ++j)
          for (size_type i = j; i < M_; ++i) {
            adjL.coeffRef(i, j) = _variRefL[++pos]->adj_;
            LA.coeffRef(i, j) = L_[pos];
          }

        for (int i = M_ - 1; i >= 0; --i) {
          for (int j = i; j >= 0; --j) {
            if (i == j) 
              adjA.coeffRef(i, j) = 0.5 * adjL.coeffRef(i, j) / LA.coeffRef(i, j);
            else {
              adjA.coeffRef(i, j) = adjL.coeffRef(i, j) / LA.coeffRef(j, j);
              adjL.coeffRef(j, j) = adjL.coeffRef(j, j) - adjL.coeffRef(i, j) * LA.coeffRef(i, j) / LA.coeffRef(j, j);
            }
            for (int k = j - 1; k >=0; --k) {
              adjL.coeffRef(i, k) = adjL.coeffRef(i, k) - adjA.coeffRef(i, j) * LA.coeffRef(j, k);
              adjL.coeffRef(j, k) = adjL.coeffRef(j, k) - adjA.coeffRef(i, j) * LA.coeffRef(i, k);
            }
          }
        }
        pos = 0;
        for (size_type j = 0; j < M_; ++j)
          for (size_type i = j; i < M_; ++i) 
            _variRefA[pos++]->adj_ += adjA.coeffRef(i, j);
      }
    };

    Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic>
    cholesky_decompose(const Eigen::Matrix<var, -1, -1> &A) {
      stan::math::check_square("cholesky_decompose", "A", A);
      stan::math::check_symmetric("cholesky_decompose", "A", A);
      
      // NOTE: this is not a memory leak, this vari is used in the
      // expression graph to evaluate the adjoint, but is not needed
      // for the returned matrix.  Memory will be cleaned up with the
      // arena allocator.
      cholesky_decompose_v_vari *baseVari
        = new cholesky_decompose_v_vari(A);

      Eigen::Matrix<var, -1, -1> L(A.rows(),A.rows());
      size_t pos = 0;
      for (size_type j = 0; j < L.cols(); ++j) {
        for (size_type i = j; i < L.rows(); ++i) 
          L.coeffRef(i, j).vi_ = baseVari->_variRefL[pos++];
        for (size_type k = 0; k < j; ++k) 
          L.coeffRef(k, j).vi_ = baseVari->_dummy[0];
      }
      return L;
    }
  }
}
#endif
