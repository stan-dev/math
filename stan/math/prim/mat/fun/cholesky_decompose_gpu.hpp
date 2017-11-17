#ifndef STAN_MATH_PRIM_MAT_FUN_CHOLESKY_DECOMPOSE_GPU_HPP
#define STAN_MATH_PRIM_MAT_FUN_CHOLESKY_DECOMPOSE_GPU_HPP

#include <stan/math/prim/mat/fun/ocl_gpu.hpp>
#include <stan/math/prim/arr/fun/matrix_gpu.hpp>
#include <stan/math/prim/mat/fun/basic_matrix_gpu.hpp>
#include <stan/math/prim/mat/fun/multiply_gpu.hpp>
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
      check_symmetric("cholesky_decompose", "m", m);
      Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
       m_tmp(m.rows(), m.cols());
      check_symmetric("cholesky_decompose", "m", m);
      matrix_gpu A(m);
      cl::Kernel kernel_chol_block = get_kernel("cholesky_block");
      cl::Kernel kernel_zero = get_kernel("cholesky_zero");
      cl::CommandQueue cmd_queue = get_queue();
      try {
        cl::Context& ctx = get_context();
        // Will be managed by the library core system
        int block = 64;
        int offset = 0;
        matrix_gpu V(block, block);
        cl::Buffer buffer_D(ctx, CL_MEM_READ_WRITE,
         sizeof(T) * block * block);
        kernel_chol_block.setArg(0, A.buffer());
        kernel_chol_block.setArg(1, offset);
        kernel_chol_block.setArg(2, A.rows());
        kernel_chol_block.setArg(3, block);
        kernel_chol_block.setArg(5, buffer_D);

        kernel_zero.setArg(0, A.buffer());
        kernel_zero.setArg(1, A.rows());

        matrix_gpu Mid;
        while ((offset + block) < (A.rows())) {
          matrix_gpu L(A.rows()-offset-block, block);
          matrix_gpu L1(A.rows()-offset-block, block);
          matrix_gpu Mid_temp(A.rows()-offset-block, A.rows()-offset-block);

          kernel_chol_block.setArg(1, offset);
          kernel_chol_block.setArg(4, V.buffer());
          cmd_queue.enqueueNDRangeKernel(kernel_chol_block,
           cl::NullRange, cl::NDRange(block), cl::NDRange(block));

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
          kernel_chol_block.setArg(4, V.buffer());
          kernel_chol_block.setArg(1, offset);
          kernel_chol_block.setArg(3, left);
          cmd_queue.enqueueNDRangeKernel(kernel_chol_block,
           cl::NullRange, cl::NDRange(left), cl::NDRange(left));
        }
        cmd_queue.enqueueNDRangeKernel(kernel_zero,
         cl::NullRange, cl::NDRange(A.rows(), A.rows()), cl::NullRange);
      } catch (const cl::Error& e) {
        check_ocl_error("cholesky_decompose", e);
      }
      copy_triangular_transposed(A, LOWER_TO_UPPER_TRIANGULAR);
      check_positive_definite_gpu("cholesky_decompose_gpu",
        "Matrix m", A);
      copy(A, m_tmp); // NOLINT
      return m_tmp.template triangularView<Eigen::Lower>();
    }
  }
}

#endif
