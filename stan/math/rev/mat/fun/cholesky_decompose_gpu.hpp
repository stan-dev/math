#ifndef STAN_MATH_REV_MAT_FUN_CHOLESKY_DECOMPOSE_GPU_HPP
#define STAN_MATH_REV_MAT_FUN_CHOLESKY_DECOMPOSE_GPU_HPP

#ifdef STAN_GPU
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/arr/fun/matrix_gpu.hpp>
#include <stan/math/prim/mat/fun/basic_matrix_gpu.hpp>
#include <stan/math/prim/mat/fun/cholesky_decompose_gpu.hpp>
#include <stan/math/prim/mat/fun/multiply_gpu.hpp>
#include <stan/math/prim/mat/fun/inverse_gpu.hpp>
#include <stan/math/prim/mat/fun/typedefs.hpp>
#include <stan/math/rev/scal/fun/value_of_rec.hpp>
#include <stan/math/rev/scal/fun/value_of.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/mat/fun/value_of_rec.hpp>
#include <stan/math/prim/mat/err/check_pos_definite.hpp>
#include <stan/math/prim/mat/err/check_square.hpp>
#include <stan/math/prim/mat/err/check_symmetric.hpp>
#include <algorithm>
#include <vector>

namespace stan {
  namespace math {

    class cholesky_gpu : public vari {
    public:
      int M_;
      vari** variRefA_;
      vari** variRefL_;


      /**
       * Constructor for GPU cholesky function.
       *
       * Stores varis for A.  Instantiates and stores varis for L.
       * Instantiates and stores dummy vari for upper triangular part of var
       * result returned in cholesky_decompose function call
       *
       * variRefL aren't on the chainable autodiff stack, only used for storage
       * and computation. Note that varis for L are constructed externally in
       * cholesky_decompose.
       *
       *
       * @param A matrix
       * @param L_A matrix, cholesky factor of A
       */
      cholesky_gpu(const Eigen::Matrix<var, -1, -1>& A,
                   const Eigen::Matrix<double, -1, -1>& L_A)
        : vari(0.0),
          M_(A.rows()),
          variRefA_(ChainableStack::memalloc_.alloc_array<vari*>
                    (A.rows() * (A.rows() + 1) / 2)),
          variRefL_(ChainableStack::memalloc_.alloc_array<vari*>
                    (A.rows() * (A.rows() + 1) / 2)) {
            size_t pos = 0;
            for (size_type j = 0; j < M_; ++j) {
              for (size_type i = j; i < M_; ++i) {
                variRefA_[pos] = A.coeffRef(i, j).vi_;
                variRefL_[pos] = new vari(L_A.coeffRef(i, j), false); ++pos;
              }
            }
          }

      /**
       * Reverse mode differentiation algorithm using a GPU
       * 
       * Reference:
       *
       * Iain Murray: Differentiation of the Cholesky decomposition, 2016.
       *
       */
      virtual void chain() {
        using Eigen::MatrixXd;
        using Eigen::Lower;
        using Eigen::Block;
        using Eigen::Upper;
        using Eigen::StrictlyUpper;
        using Eigen::StrictlyLower;
        MatrixXd Lbar(M_, M_);
        MatrixXd Lbar1(M_, M_);
        MatrixXd L(M_, M_);
        Lbar.setZero();
        L.setZero();

        size_t pos = 0;
        for (size_type j = 0; j < M_; ++j) {
          for (size_type i = j; i < M_; ++i) {
            Lbar.coeffRef(i, j) = variRefL_[pos]->adj_;
            L.coeffRef(i, j) = variRefL_[pos]->val_;
            ++pos;
          }
        }
/*
        write_binary("~/L10000.dat",L);
        write_binary("~/Lbar10000.dat",Lbar);
        exit(1);*/
        matrix_gpu L_gpu(L);
        matrix_gpu Lbar_gpu(Lbar);

        int M = M_;
        int block_size_ = 128;
        block_size_ = std::max((M / 8 / 16) * 16, 8);
        block_size_ = std::min(block_size_, 512);
        for (int k = M; k > 0; k -= block_size_) {
          int j = std::max(0, k - block_size_);
          matrix_gpu R_gpu(k-j, j);
          matrix_gpu D_gpu(k-j, k-j);
          matrix_gpu Dinv_gpu(k-j, k-j);
          matrix_gpu B_gpu(M-k, j);
          matrix_gpu C_gpu(M-k, k-j);

          matrix_gpu Rbar_gpu(k-j, j);
          matrix_gpu Dbar_gpu(k-j, k-j);
          matrix_gpu Dbar2_gpu(k-j, k-j);
          matrix_gpu Bbar_gpu(M-k, j);
          matrix_gpu Bbar2_gpu(M-k, j);
          matrix_gpu Cbar_gpu(M-k, k-j);
          matrix_gpu Cbar2_gpu(k-j, M-k);
          matrix_gpu Cbar3_gpu(k-j, M-k);
          matrix_gpu temp_gpu(k-j, j);

          copy_submatrix(L_gpu, R_gpu, j, 0, 0, 0, k-j, j);
          copy_submatrix(L_gpu, D_gpu, j, j, 0, 0, k-j, k-j);
          copy_submatrix(L_gpu, B_gpu, k, 0, 0, 0, M-k, j);
          copy_submatrix(L_gpu, C_gpu, k, j, 0, 0, M-k, k-j);

          copy_submatrix(Lbar_gpu, Rbar_gpu, j, 0, 0, 0, k-j, j);
          copy_submatrix(Lbar_gpu, Dbar_gpu, j, j, 0, 0, k-j, k-j);
          copy_submatrix(Lbar_gpu, Bbar_gpu, k, 0, 0, 0, M-k, j);
          copy_submatrix(Lbar_gpu, Cbar_gpu, k, j, 0, 0, M-k, k-j);

          if (Cbar_gpu.size() > 0) {
            stan::math::copy(D_gpu, Dinv_gpu);
            lower_triangular_inverse(Dinv_gpu);
            Dinv_gpu = transpose(Dinv_gpu);
            Cbar_gpu = transpose(Cbar2_gpu);

            Cbar3_gpu = multiply(Dinv_gpu, Cbar2_gpu);
            Cbar3_gpu = transpose(Cbar_gpu);

            Bbar2_gpu = multiply(Cbar_gpu, R_gpu);
            Bbar_gpu = subtract(Bbar_gpu, Bbar2_gpu);

            Cbar_gpu = transpose(Cbar3_gpu);
            Dbar2_gpu = multiply(Cbar3_gpu, C_gpu);
            Dbar_gpu = subtract(Dbar_gpu, Dbar2_gpu);
          }

          D_gpu = transpose(D_gpu);
          zeros(Dbar_gpu, UPPER);
          Dbar2_gpu = multiply(D_gpu, Dbar_gpu);
          copy_triangular_transposed(Dbar2_gpu,
            LOWER_TO_UPPER_TRIANGULAR);
          D_gpu = transpose(D_gpu);
          lower_triangular_inverse(D_gpu);
          D_gpu = transpose(D_gpu);
          Dbar_gpu = multiply(D_gpu, Dbar2_gpu);
          Dbar_gpu = transpose(Dbar_gpu);
          Dbar2_gpu = multiply(D_gpu, Dbar_gpu);

          if (Cbar_gpu.size() > 0 && B_gpu.size() > 0) {
            Cbar_gpu = transpose(Cbar2_gpu);
            temp_gpu = multiply(Cbar2_gpu, B_gpu);
            Rbar_gpu = subtract(Rbar_gpu, temp_gpu);
          }

          if (Dbar_gpu.size() > 0 && R_gpu.size() > 0) {
            stan::math::copy(Dbar2_gpu, Dbar_gpu);
            copy_triangular_transposed(Dbar_gpu,
              LOWER_TO_UPPER_TRIANGULAR);
            temp_gpu = multiply(Dbar_gpu, R_gpu);
            Rbar_gpu = subtract(Rbar_gpu, temp_gpu);
          }

          diagonal_multiply(Dbar2_gpu, 0.5);
          zeros(Dbar2_gpu, UPPER);

          copy_submatrix(Rbar_gpu, Lbar_gpu, 0, 0, j, 0, k-j, j);
          copy_submatrix(Dbar2_gpu, Lbar_gpu, 0, 0, j, j, k-j, k-j);
          copy_submatrix(Bbar_gpu, Lbar_gpu, 0, 0, k, 0, M-k, j);
          copy_submatrix(Cbar_gpu, Lbar_gpu, 0, 0, k, j, M-k, k-j);
        }
        stan::math::copy(Lbar_gpu, Lbar);
        pos = 0;
        for (size_type j = 0; j < M_; ++j)
          for (size_type i = j; i < M_; ++i)
            variRefA_[pos++]->adj_ += Lbar.coeffRef(i, j);
      }
    };
    /**
     * Reverse mode specialization of cholesky decomposition on a GPU
     *
     * Note chainable stack varis are created below in Matrix<var, -1, -1>
     *
     * @param A Matrix
     * @return L cholesky factor of A
     */
    inline Eigen::Matrix<var, -1, -1>
      cholesky_decompose_gpu(const Eigen::Matrix<var, -1, -1> &A) {
      check_square("cholesky_decompose", "A", A);
      check_symmetric("cholesky_decompose", "A", A);

      Eigen::Matrix<double, -1, -1> L_A(value_of_rec(A));
      L_A = L_A.selfadjointView<Eigen::Lower>();
      L_A = cholesky_decompose_gpu(L_A);
      // Memory allocated in arena.
      vari* dummy = new vari(0.0, false);
      Eigen::Matrix<var, -1, -1> L(A.rows(), A.cols());
      cholesky_gpu *baseVari = new cholesky_gpu(A, L_A);
      size_t pos = 0;
      for (size_type j = 0; j < L.cols(); ++j) {
        for (size_type i = j; i < L.cols(); ++i) {
          L.coeffRef(i, j).vi_ = baseVari->variRefL_[pos++];
        }
        for (size_type k = 0; k < j; ++k)
          L.coeffRef(k, j).vi_ = dummy;
      }
      return L;
    }
  }
}
#endif
#endif
