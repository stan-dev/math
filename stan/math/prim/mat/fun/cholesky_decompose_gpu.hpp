#ifndef STAN_MATH_PRIM_MAT_FUN_CHOLESKY_DECOMPOSE_GPU_HPP
#define STAN_MATH_PRIM_MAT_FUN_CHOLESKY_DECOMPOSE_GPU_HPP

#include <stan/math/prim/mat/fun/ocl_gpu.hpp>
#include <stan/math/prim/arr/fun/matrix_gpu.hpp>
#include <stan/math/prim/mat/fun/basic_matrix_gpu.hpp>
#include <stan/math/prim/mat/fun/multiply_gpu.hpp>
#include <stan/math/prim/mat/fun/inverse_gpu.hpp>
#include <stan/math/prim/mat/err/check_gpu.hpp>
#include <stan/math/prim/mat/err/check_symmetric.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/rev/scal/fun/value_of_rec.hpp>
#include <stan/math/rev/scal/fun/value_of.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_arithmetic.hpp>
#include <iostream>
#include <string>
#include <map>

/*   @file stanmathcl/matrix_inverse.hpp
*    @brief matrix_inverse  -  functions for matrix inversion:
*     lower triangular,  upper triangular,  regular,  ...
*/

// CURRENTLY ONLY SUPPORTS LOWER TRIANGULAR
namespace stan {
  namespace math {
    /**
     * Return the lower-triangular Cholesky factor (i.e., matrix
     * square root) of the specified square, symmetric matrix.  
     * The return value \f$L\f$ will be a lower-traingular matrix such that the
     * original matrix \f$A\f$ is given by
     * <p>\f$A = L \times L^T\f$.
     * The Cholesky decomposition is computed on the GPU. The
     * input matrix is transfered to the GPU and the resulting
     * lower-triangular matrix is then copied from the GPU.
     * 
     * @param A Symmetrix matrix on the GPU.
     * @param block size of the block for each step
     * @return Square root of matrix on the GPU.
     * @throw std::domain_error if m is not
     *  positive definite (if m has more than 0 elements)
     */    
    inline matrix_gpu cholesky_decompose_gpu(matrix_gpu& A, int block) {
      cl::Kernel kernel_chol_block = get_kernel("cholesky_block");
      cl::CommandQueue cmd_queue = get_queue();

      // Will be managed by the library core system
      int offset = 0;
      matrix_gpu V(block, block);
      matrix_gpu D(block, block);
      while ((offset + block) < (A.rows())) {
        matrix_gpu L(A.rows()-offset-block, block);
        matrix_gpu Mid(A.rows()-offset-block, A.rows()-offset-block);
        matrix_gpu Mid_temp(A.rows()-offset-block, A.rows()-offset-block);

        copy_submatrix(A, D, offset, offset, 0, 0, block, block);
        zeros(V);
        try {
          kernel_chol_block.setArg(0, V.buffer());
          kernel_chol_block.setArg(1, D.buffer());
          kernel_chol_block.setArg(2, block);
          cmd_queue.enqueueNDRangeKernel(kernel_chol_block,
           cl::NullRange, cl::NDRange(block), cl::NDRange(block));
        } catch (const cl::Error& e) {
          check_ocl_error("cholesky_decompose", e);
        }
        copy(V, D);
        copy_submatrix(V, A, 0, 0, offset, offset, block, block);

        V = lower_triangular_inverse(D);

        copy_submatrix(A, L, (offset+block), offset, 0, 0,
          (A.rows()-offset-block) , block);
        V = transpose(V);
        L = multiply(L, V);
        copy_submatrix(L, A, 0, 0, (offset+block), offset,
          (A.rows()-offset-block) , block);

        copy_submatrix(A, Mid_temp, (offset+block), (offset+block),
          0, 0, (A.rows()-offset-block), (A.rows()-offset-block));
        Mid = multiply_with_self_transposed(L);
        Mid = subtract(Mid_temp, Mid);
        copy_submatrix(Mid, A, 0, 0, (offset+block), (offset+block),
          (A.rows()-offset-block), (A.rows()-offset-block));
        offset += block;
      }
      int left = A.rows() - offset;
      if (left > 0) {
        matrix_gpu D(left, left);
        matrix_gpu V(left, left);
        copy_submatrix(A, D, offset, offset, 0, 0, left, left);
        zeros(V);
        try {
          kernel_chol_block.setArg(0, V.buffer());
          kernel_chol_block.setArg(1, D.buffer());
          kernel_chol_block.setArg(2, left);
          cmd_queue.enqueueNDRangeKernel(kernel_chol_block,
           cl::NullRange, cl::NDRange(left), cl::NDRange(left));
        } catch (const cl::Error& e) {
          check_ocl_error("cholesky_decompose", e);
        }
        copy_submatrix(V, A, 0, 0, offset, offset, left, left);
      }
      zeros(A, UPPER);
      copy_triangular_transposed(A, LOWER_TO_UPPER_TRIANGULAR);
      check_positive_definite("cholesky_decompose_gpu",
        "Matrix m", A);
      zeros(A, UPPER);
      matrix_gpu B(A);
      return B;
    }
     /**
     * Return the lower-triangular Cholesky factor (i.e., matrix
     * square root) of the specified square, symmetric matrix.  
     * The return value \f$L\f$ will be a lower-traingular matrix such that the
     * original matrix \f$A\f$ is given by
     * <p>\f$A = L \times L^T\f$.
     * The Cholesky decomposition is computed on the GPU. The
     * input matrix is transfered to the GPU and the resulting
     * lower-triangular matrix is then copied from the GPU.
     * 
     * @param m Symmetrix matrix.
     * @return Square root of matrix.
     * @throw std::domain_error if m is not a symmetric matrix or
     *   if m is not positive definite (if m has more than 0 elements)
     */
    template <typename T>
    typename boost::enable_if_c<boost::is_arithmetic<T>::value,
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>::type
    cholesky_decompose_gpu(const Eigen::Matrix<T,
     Eigen::Dynamic, Eigen::Dynamic>& m) {
            if (m.size() == 0) return m;

      matrix_gpu A(m);

      check_symmetric("cholesky_decompose", "m", A);

      Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
       m_tmp(m.rows(), m.cols());

      cl::CommandQueue cmd_queue = get_queue();
      cl::Kernel kernel_chol_block = get_kernel("cholesky_block");

      // Will be managed by the library core system
      int block = 400;
      int offset = 0;
      matrix_gpu V(block, block);
      matrix_gpu D(block, block);

      while ((offset + block) < (A.rows())) {
        matrix_gpu L(A.rows()-offset-block, block);
        matrix_gpu Mid(A.rows()-offset-block, A.rows()-offset-block);
        matrix_gpu Mid_temp(A.rows()-offset-block, A.rows()-offset-block);

        copy_submatrix(A, D, offset, offset, 0, 0, block, block);
        zeros(V);
        V = cholesky_decompose_gpu(D, 100);
        copy(V, D);
        copy_submatrix(V, A, 0, 0, offset, offset, block, block);

        V = lower_triangular_inverse(D);

        copy_submatrix(A, L, (offset+block), offset, 0, 0,
          (A.rows()-offset-block) , block);
        V = transpose(V);
        L = multiply(L, V);
        copy_submatrix(L, A, 0, 0, (offset+block), offset,
          (A.rows()-offset-block) , block);

        copy_submatrix(A, Mid_temp, (offset+block), (offset+block),
          0, 0, (A.rows()-offset-block), (A.rows()-offset-block));
        Mid = multiply_with_self_transposed(L);
        Mid = subtract(Mid_temp, Mid);
        copy_submatrix(Mid, A, 0, 0, (offset+block), (offset+block),
          (A.rows()-offset-block), (A.rows()-offset-block));

        offset += block;
      }
      int left = A.rows() - offset;
      if (left > 0) {
        matrix_gpu D(left, left);
        matrix_gpu V(left, left);
        copy_submatrix(A, D, offset, offset, 0, 0, left, left);
        zeros(V);
        V = cholesky_decompose_gpu(D, 100);
        copy_submatrix(V, A, 0, 0, offset, offset, left, left);
      }
      zeros(A, UPPER);
      copy_triangular_transposed(A, LOWER_TO_UPPER_TRIANGULAR);
      check_positive_definite("cholesky_decompose_gpu",
        "Matrix m", A);
      zeros(A, UPPER);
      copy(A, m_tmp); // NOLINT
      return m_tmp;
    }
  }
}

#endif
