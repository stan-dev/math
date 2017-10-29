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
          stan::math::matrix_gpu R_gpu(k-j, j);
          stan::math::matrix_gpu D_gpu(k-j, k-j);
          stan::math::matrix_gpu Dinv_gpu(k-j, k-j);
          stan::math::matrix_gpu B_gpu(M-k, j);
          stan::math::matrix_gpu C_gpu(M-k, k-j);

          stan::math::matrix_gpu Rbar_gpu(k-j, j);
          stan::math::matrix_gpu Dbar_gpu(k-j, k-j);
          stan::math::matrix_gpu Dbar2_gpu(k-j, k-j);
          stan::math::matrix_gpu Bbar_gpu(M-k, j);
          stan::math::matrix_gpu Bbar2_gpu(M-k, j);
          stan::math::matrix_gpu Cbar_gpu(M-k, k-j);
          stan::math::matrix_gpu Cbar2_gpu(k-j, M-k);
          stan::math::matrix_gpu Cbar3_gpu(k-j, M-k);
          stan::math::matrix_gpu temp_gpu(k-j, j);

          stan::math::copy_submatrix(L_gpu, R_gpu, j, 0, 0, 0, k-j, j);
          stan::math::copy_submatrix(L_gpu, D_gpu, j, j, 0, 0, k-j, k-j);
          stan::math::copy_submatrix(L_gpu, B_gpu, k, 0, 0, 0, M-k, j);
          stan::math::copy_submatrix(L_gpu, C_gpu, k, j, 0, 0, M-k, k-j);

          stan::math::copy_submatrix(Lbar_gpu, Rbar_gpu, j, 0, 0, 0, k-j, j);
          stan::math::copy_submatrix(Lbar_gpu, Dbar_gpu, j, j, 0, 0, k-j, k-j);
          stan::math::copy_submatrix(Lbar_gpu, Bbar_gpu, k, 0, 0, 0, M-k, j);
          stan::math::copy_submatrix(Lbar_gpu, Cbar_gpu, k, j, 0, 0, M-k, k-j);

          if (Cbar_gpu.size() > 0) {
            stan::math::copy(D_gpu, Dinv_gpu);
            stan::math::lower_triangular_inverse(Dinv_gpu);
            stan::math::transpose(Dinv_gpu);
            stan::math::transpose(Cbar2_gpu, Cbar_gpu);

            stan::math::multiply(Dinv_gpu, Cbar2_gpu, Cbar3_gpu);
            stan::math::transpose(Cbar_gpu, Cbar3_gpu);

            stan::math::multiply(Cbar_gpu, R_gpu, Bbar2_gpu);
            stan::math::subtract(Bbar_gpu, Bbar_gpu, Bbar2_gpu);

            stan::math::transpose(Cbar3_gpu, Cbar_gpu);
            stan::math::multiply(Cbar3_gpu, C_gpu, Dbar2_gpu);
            stan::math::subtract(Dbar_gpu, Dbar_gpu, Dbar2_gpu);
          }

          stan::math::transpose(D_gpu);
          stan::math::zeros(Dbar_gpu, stan::math::UPPER);
          stan::math::multiply(D_gpu, Dbar_gpu, Dbar2_gpu);
          stan::math::copy_triangular_transposed(Dbar2_gpu,
            stan::math::LOWER_TO_UPPER_TRIANGULAR);
          stan::math::transpose(D_gpu);
          stan::math::lower_triangular_inverse(D_gpu);
          stan::math::transpose(D_gpu);
          stan::math::multiply(D_gpu, Dbar2_gpu, Dbar_gpu);
          stan::math::transpose(Dbar_gpu);
          stan::math::multiply(D_gpu, Dbar_gpu, Dbar2_gpu);

          if (Cbar_gpu.size() > 0 && B_gpu.size() > 0) {
            stan::math::transpose(Cbar2_gpu, Cbar_gpu);
            stan::math::multiply(Cbar2_gpu, B_gpu, temp_gpu);
            stan::math::subtract(Rbar_gpu, Rbar_gpu, temp_gpu);
          }

          if (Dbar_gpu.size() > 0 && R_gpu.size() > 0) {
            stan::math::copy(Dbar2_gpu, Dbar_gpu);
            stan::math::copy_triangular_transposed(Dbar_gpu,
              stan::math::LOWER_TO_UPPER_TRIANGULAR);
            stan::math::multiply(Dbar_gpu, R_gpu, temp_gpu);
            stan::math::subtract(Rbar_gpu, Rbar_gpu, temp_gpu);
          }

          stan::math::diagonal_multiply(Dbar2_gpu, 0.5);
          stan::math::zeros(Dbar2_gpu, stan::math::UPPER);

          stan::math::copy_submatrix(Rbar_gpu, Lbar_gpu, 0, 0, j, 0, k-j, j);
          stan::math::copy_submatrix(Dbar2_gpu, Lbar_gpu, 0, 0, j, j, k-j, k-j);
          stan::math::copy_submatrix(Bbar_gpu, Lbar_gpu, 0, 0, k, 0, M-k, j);
          stan::math::copy_submatrix(Cbar_gpu, Lbar_gpu, 0, 0, k, j, M-k, k-j);
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
      L_A = stan::math::cholesky_decompose_gpu(L_A);
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
