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
       * Reverse mode differentiation algorithm refernce:
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
        stan::math::matrix_gpu L_gpu(L);
        stan::math::matrix_gpu Lbar_gpu(Lbar);  
        stan::math::matrix_gpu Lbar_temp_gpu(M_,M_);
        
        /***The CPU version of the derivative **/
        /* If this is here and the result of this CPU derivative is not used, the tests pass, otherwise NaN is returend for some weird reason*/
        /* */        
        L.transposeInPlace();
        Lbar = (L * Lbar);
        
        Lbar.triangularView<Eigen::StrictlyUpper>() = Lbar.transpose().triangularView<Eigen::StrictlyUpper>();
        L.triangularView<Eigen::Upper>().solveInPlace(Lbar);
        L.triangularView<Eigen::Upper>().solveInPlace(Lbar.transpose());
        
        Lbar.transposeInPlace();
        Lbar.diagonal() *= 0.5;
        /****/
        
        stan::math::transpose(L_gpu);
        stan::math::multiply(L_gpu, Lbar_gpu, Lbar_temp_gpu);
        stan::math::copy_triangular_transposed(Lbar_temp_gpu,
              LOWER_TO_UPPER_TRIANGULAR);        
        stan::math::transpose(L_gpu);
        stan::math::lower_triangular_inverse(L_gpu);
        stan::math::transpose(L_gpu);
        stan::math::transpose(Lbar_temp_gpu);
        stan::math::multiply(L_gpu,Lbar_temp_gpu,Lbar_gpu);
        stan::math::transpose(Lbar_gpu);
        stan::math::multiply(L_gpu,Lbar_gpu,Lbar_temp_gpu);
        stan::math::transpose(Lbar_temp_gpu);
        stan::math::diagonal_multiply_with_scalar(Lbar_temp_gpu,0.5);
        stan::math::copy(Lbar_temp_gpu, Lbar1);
        /*
        double err=0, diff=0;
        for(int i=0;i<M_;i++){
          for(int j=0;j<M_;j++){
            diff=fabs(Lbar1(i,j)-Lbar(i,j));
            if(diff>err){
              err=diff;
            }
          }
        }
        std::cout << "MERR: " << err << std::endl;
        */
        pos = 0;
        for (size_type j = 0; j < M_; ++j)
          for (size_type i = j; i < M_; ++i)
            variRefA_[pos++]->adj_ += Lbar1.coeffRef(i, j);
      }
    };
    /**
     * Reverse mode specialization of cholesky decomposition
     *
     * Internally calls Eigen::LLT rather than using
     * stan::math::cholesky_decompose in order to use an inplace decomposition.
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
//      std::cout << "OUTPUT REV: \n";
//    TODO: Have the matrix stay on the GPU for this and the derivative
      L_A = stan::math::cholesky_decompose_gpu(L_A);
//      std::cout << "OUTPUT REV: \n";
//      std::cout << L_A << "\n";
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
