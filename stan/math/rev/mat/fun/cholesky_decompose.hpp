#ifndef STAN_MATH_REV_MAT_FUN_CHOLESKY_DECOMPOSE_HPP
#define STAN_MATH_REV_MAT_FUN_CHOLESKY_DECOMPOSE_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/fun/typedefs.hpp>
#include <stan/math/prim/mat/fun/cholesky_decompose.hpp>
#include <stan/math/rev/scal/fun/value_of_rec.hpp>
#include <stan/math/rev/scal/fun/value_of.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/mat/fun/value_of_rec.hpp>
#include <stan/math/prim/mat/err/check_pos_definite.hpp>
#include <stan/math/prim/mat/err/check_square.hpp>
#include <stan/math/prim/mat/err/check_symmetric.hpp>

namespace stan {
  namespace math {

    class cholesky_decompose_v_vari : public vari {
    public:
      int M_;  // A.rows() = A.cols()
      vari** variRefA_;
      vari** variRefL_;

      /* ctor for cholesky function
       *
       * Stores varis for A
       * Instantiates and stores varis for L
       * Instantiates and stores dummy vari for
       * upper triangular part of var result returned
       * in cholesky_decompose function call
       *
       * variRefL aren't on the chainable
       * autodiff stack, only used for storage
       * and computation. Note that varis for
       * L are constructed externally in
       * cholesky_decompose.
       *
       * @param matrix A
       * @param matrix L, cholesky factor of A
       * */
      cholesky_decompose_v_vari(const Eigen::Matrix<var, -1, -1>& A,
                                const Eigen::Matrix<double, -1, -1>& L_A)
        : vari(0.0),
          M_(A.rows()),
          variRefA_(ChainableStack::memalloc_.alloc_array<vari*>
                    (A.rows() * (A.rows() + 1) / 2)),
          variRefL_(ChainableStack::memalloc_.alloc_array<vari*>
                    (A.rows() * (A.rows() + 1) / 2)) {
        size_t accum = 0;
        size_t accum_i = accum;
        for (size_type j = 0; j < M_; ++j) {
          for (size_type i = j; i < M_; ++i) {
            accum_i += i;
            size_t pos = j + accum_i;
            variRefA_[pos] = A.coeffRef(i, j).vi_;
            variRefL_[pos] = new vari(L_A.coeffRef(i, j), false);
          }
          accum += j;
          accum_i = accum;
        }
      }

      /* Reverse mode differentiation
       * algorithm refernce:
       *
       * Mike Giles. An extended collection of matrix
       * derivative results for forward and reverse mode AD.
       * Jan. 2008.
       *
       * Note algorithm  as laid out in Giles is
       * row-major, so Eigen::Matrices are explicitly storage
       * order RowMajor, whereas Eigen defaults to
       * ColumnMajor. Also note algorithm
       * starts by calculating the adjoint for
       * A(M_ - 1, M_ - 1), hence pos on line 94 is decremented
       * to start at pos = M_ * (M_ + 1) / 2.
       * */
      virtual void chain() {
        using Eigen::Matrix;
        using Eigen::RowMajor;
        Matrix<double, -1, -1, RowMajor> adjL(M_, M_);
        Matrix<double, -1, -1, RowMajor> LA(M_, M_);
        Matrix<double, -1, -1, RowMajor> adjA(M_, M_);
        size_t pos = 0;
        for (size_type i = 0; i < M_; ++i) {
          for (size_type j = 0; j <= i; ++j) {
            adjL.coeffRef(i, j) = variRefL_[pos]->adj_;
            LA.coeffRef(i, j) = variRefL_[pos]->val_;
            ++pos;
          }
        }

        --pos;
        for (int i = M_ - 1; i >= 0; --i) {
          for (int j = i; j >= 0; --j) {
            if (i == j) {
              adjA.coeffRef(i, j) = 0.5 * adjL.coeff(i, j)
                / LA.coeff(i, j);
            } else {
              adjA.coeffRef(i, j) = adjL.coeff(i, j)
                / LA.coeff(j, j);
              adjL.coeffRef(j, j) -= adjL.coeff(i, j)
                * LA.coeff(i, j) / LA.coeff(j, j);
            }
            for (int k = j - 1; k >=0; --k) {
              adjL.coeffRef(i, k) -= adjA.coeff(i, j)
                * LA.coeff(j, k);
              adjL.coeffRef(j, k) -= adjA.coeff(i, j)
                * LA.coeff(i, k);
            }
            variRefA_[pos--]->adj_ += adjA.coeffRef(i, j);
          }
        }
      }
    };

    /* Reverse mode specialization of
     * cholesky decomposition
     *
     * Internally calls llt rather than using
     * cholesky_decompose in order
     * to use selfadjointView<Lower> optimization.
     *
     * Note chainable stack varis are created
     * below in Matrix<var, -1, -1>
     *
     * @param Matrix A
     * @return L cholesky factor of A
     */
    inline Eigen::Matrix<var, -1, -1>
      cholesky_decompose(const Eigen::Matrix<var, -1, -1> &A) {
      check_square("cholesky_decompose", "A", A);
      check_symmetric("cholesky_decompose", "A", A);

      Eigen::Matrix<double, -1, -1> L_A(value_of_rec(A));
      Eigen::LLT<Eigen::MatrixXd> L_factor
        = L_A.selfadjointView<Eigen::Lower>().llt();
      check_pos_definite("cholesky_decompose", "m", L_factor);
      L_A = L_factor.matrixL();

      // NOTE: this is not a memory leak, this vari is used in the
      // expression graph to evaluate the adjoint, but is not needed
      // for the returned matrix.  Memory will be cleaned up with the
      // arena allocator.
      cholesky_decompose_v_vari *baseVari
        = new cholesky_decompose_v_vari(A, L_A);
      vari* dummy = new vari(0.0, false);
      Eigen::Matrix<var, -1, -1> L(A.rows(), A.cols());
      size_t accum = 0;
      size_t accum_i = accum;
      for (size_type j = 0; j < L.cols(); ++j) {
        for (size_type i = j; i < L.cols(); ++i) {
          accum_i += i;
          size_t pos = j + accum_i;
          L.coeffRef(i, j).vi_ = baseVari->variRefL_[pos];
        }
        for (size_type k = 0; k < j; ++k)
          L.coeffRef(k, j).vi_ = dummy;
        accum += j;
        accum_i = accum;
      }
      return L;
    }

  }
}
#endif
