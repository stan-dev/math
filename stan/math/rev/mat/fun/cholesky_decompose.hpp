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

    class cholesky_decompose_v_vari : public vari {
    public:
      int M_;  // A.rows() = A.cols() = B.rows()
      double* L_;
      vari** _variRefA;
      vari** _variRefL;
      vari** _dummy;

      cholesky_decompose_v_vari(const Eigen::Matrix<var, -1, -1>& A,
                                const Eigen::Matrix<double, -1, -1>& L_A)
        : vari(0.0),
          M_(A.rows()),
          L_(reinterpret_cast<double*>
             (stan::math::ChainableStack::memalloc_
              .alloc(sizeof(double) * A.rows() * (A.rows() + 1) / 2))),
          _variRefA(reinterpret_cast<vari**>
                    (stan::math::ChainableStack::memalloc_
                     .alloc(sizeof(vari*) * A.rows() * (A.cols() + 1) / 2))),
          _variRefL(reinterpret_cast<vari**>
                    (stan::math::ChainableStack::memalloc_
                     .alloc(sizeof(vari*) * A.rows() * (A.rows() + 1) / 2))),
          _dummy(reinterpret_cast<vari**>
                    (stan::math::ChainableStack::memalloc_
                     .alloc(sizeof(vari*))))  {
        using Eigen::Matrix;
        using Eigen::Map;

        size_t pos = 0;
        for (size_type j = 0; j < M_; ++j) {
          for (size_type i = j; i < M_; ++i) {
            _variRefA[pos++] = A(i, j).vi_;
          }
        }

        pos = 0;
        for (size_type j = 0; j < M_; ++j) {
          for (size_type i = j; i < M_; ++i) {
            L_[pos] = L_A(i, j);
            _variRefL[pos] = new vari(L_[pos], false);
            pos++;
          }
        }
        _dummy[0] = new vari(0.0, false);
      }

      virtual void chain() {
        using Eigen::Matrix;
        using Eigen::Dynamic;
        using Eigen::TriangularView;
        using Eigen::Lower;
        using Eigen::Map;
        Eigen::Matrix<double, Dynamic, Dynamic> adjA(M_, M_);
        Eigen::Matrix<double, Dynamic, Dynamic> adjL(M_, M_);
        Eigen::Matrix<double, Dynamic, Dynamic> LA(M_, M_);
        adjA.setZero();
        size_t pos = -1;
        for (size_type j = 0; j < M_; ++j)
          for (size_type i = j; i < M_; ++i) {
            adjL(i, j) = _variRefL[++pos]->adj_;
            LA(i, j) = L_[pos];
          }

        for (int i = M_ - 1; i >= 0; --i) {
          for (int j = i; j >= 0; --j) {
            if (i == j) 
              adjA(i, j) = 0.5 * adjL(i, j) / LA(i, j);
            else {
              adjA(i, j) = adjL(i, j) / LA(j, j);
              adjL(j, j) = adjL(j, j) - adjL(i, j) * LA(i, j) / LA(j, j);
            }
            for (int k = j - 1; k >=0; --k) {
              adjL(i, k) = adjL(i, k) - adjA(i, j) * LA(j, k);
              adjL(j, k) = adjL(j, k) - adjA(i, j) * LA(i, k);
            }
          }
        }
        pos = 0;
        for (size_type j = 0; j < M_; ++j)
          for (size_type i = j; i < M_; ++i) 
            _variRefA[pos++]->adj_ += adjA(i, j);
      }
    };

    Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic>
    cholesky_decompose(const Eigen::Matrix<var, -1, -1> &A) {
      Eigen::Matrix<var, -1, -1> L(A.rows(), A.cols());
      Eigen::Matrix<double, -1, -1> L_A(A.rows(), A.cols());
      
      L_A = cholesky_decompose(value_of_rec(A));
      // NOTE: this is not a memory leak, this vari is used in the
      // expression graph to evaluate the adjoint, but is not needed
      // for the returned matrix.  Memory will be cleaned up with the
      // arena allocator.
      cholesky_decompose_v_vari *baseVari
        = new cholesky_decompose_v_vari(A, L_A);

      size_t pos = 0;
      for (size_type j = 0; j < L.cols(); ++j) {
        for (size_type i = j; i < L.rows(); ++i) 
          L(i, j).vi_ = baseVari->_variRefL[pos++];
        for (size_type k = 0; k < j; ++k) 
          L(k, j).vi_ = baseVari->_dummy[0];
      }
      return L;
    }
  }
}
#endif
