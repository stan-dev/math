#ifndef STAN_MATH_REV_MAT_FUN_CHOLESKY_DECOMPOSE_HPP
#define STAN_MATH_REV_MAT_FUN_CHOLESKY_DECOMPOSE_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/err/check_pos_definite.hpp>
#include <stan/math/prim/mat/err/check_square.hpp>
#include <stan/math/prim/mat/err/check_symmetric.hpp>

namespace stan {
  namespace math {

    /* Index function to map row and column of matrix A to 
     * index vech(A) the half-vectorization operator (maps
     * lower triangular of symmetric matrix to vector).
     *
     * @param i row
     * @param j column
     * @param N dimension of NxN matrix
     * @param accum half-vector correction 
     * equal to sum(1:j)
     */
    inline size_t vech_indexer(int i, int j, int N, size_t accum) {
      return j * N + i - accum;
    }

    class cholesky_decompose_v_vari : public vari {
    public:
      int M_;
      double* A_;
      vari** _variRefA;
      vari** _variRefL;
      vari** _dummy;
      size_t accum_j;

      /* Ctor for cholesky function
       *
       * Performs cholesky decomposition 
       * on double version of matrix A
       * initially stored in matrix C
       *
       * variRefL aren't on the chainable
       * autodiff stack, only used for storage
       * and computation. Note that vars for
       * L are constructed externally.
       *
       * @param matrix A
       * */
      explicit cholesky_decompose_v_vari(const Eigen::Matrix<var, -1, -1>& A)
        : vari(0.0),
          M_(A.rows()),
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
                     .alloc(sizeof(vari*)))),
          accum_j(0) {
          using Eigen::Map;
          using Eigen::Matrix;

        size_t posA = 0;
        for (size_type j = 0; j < M_; ++j) {
          for (size_type k = 0; k < M_; ++k)
            A_[posA++] = A.coeffRef(k, j).vi_->val_;
          accum_j += j;
        }

        Matrix<double, -1, -1> L = Map<Matrix<double, -1, -1> > (A_, M_, M_);
        Eigen::LLT<Eigen::MatrixXd> L_factor
          = L.selfadjointView<Eigen::Lower>().llt();
        check_pos_definite("cholesky_decompose", "m", L_factor);
        L = L_factor.matrixL();

        size_t pos = 0;
        for (size_type j = 0; j < M_; ++j) {
          for (size_type i = j; i < M_; ++i) {
            _variRefA[pos] = A.coeffRef(i, j).vi_;
            _variRefL[pos] = new vari(L.coeffRef(i, j), false);
            ++pos;
          }
        }
        _dummy[0] = new vari(0.0, false);
      }

      /* Reverse mode differentiation 
       * algorithm refernce: 
       *
       * Giles. Collected matrix derivative results
       * for forward and reverse mode AD. Jan. 2008.
       *
       * sum(1:i), sum(1:j), and sum(1:k) needed for 
       * vech_indexer precomputed and decremented each
       * iteration in appropriate loop.
       * 
       * Algorithm overwrites the values in L's vari->adj_ 
       * (L's varis are stored in _variRefL)
       * to decrease speed and memory usage.
       * */
      virtual void chain() {
        size_t sum_j = accum_j;
        for (int i = M_ - 1; i >= 0; --i) {
          for (int j = i; j >= 0; --j) {
            size_t ij = vech_indexer(i, j, M_, sum_j);
            size_t jj = vech_indexer(j, j, M_, sum_j);
            if (i == j) {
             _variRefA[ij]->adj_ += 0.5
               * _variRefL[ij]->adj_ / _variRefL[ij]->val_;
            } else {
              _variRefA[ij]->adj_ += _variRefL[ij]->adj_ / _variRefL[jj]->val_;
              _variRefL[jj]->adj_ = _variRefL[jj]->adj_ - _variRefL[ij]->adj_
                * _variRefL[ij]->val_ / _variRefL[jj]->val_;
            }
            size_t sum_k = sum_j - j;
            for (int k = j - 1; k >=0; --k) {
              size_t ik = vech_indexer(i, k, M_, sum_k);
              size_t jk = vech_indexer(j, k, M_, sum_k);
              _variRefL[ik]->adj_ = _variRefL[ik]->adj_ - _variRefA[ij]->adj_
                * _variRefL[jk]->val_;
              _variRefL[jk]->adj_ = _variRefL[jk]->adj_ - _variRefA[ij]->adj_
                * _variRefL[ik]->val_;
              sum_k -= k;
            }
            sum_j -= j;
          }
          accum_j -= i;
          sum_j = accum_j;
        }
      }
    };

    /* Reverse mode specialization of
     * cholesky decomposition
     *
     * Internally calls llt rather than using 
     * stan::math::cholesky_decompose
     *
     * Note chainable stack varis are created
     * below in Matrix<var, -1, -1>.
     *
     * @param Matrix A
     */
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

      Eigen::Matrix<var, -1, -1> L(A.rows(), A.rows());
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
