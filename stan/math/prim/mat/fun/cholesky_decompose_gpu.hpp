#ifndef STAN_MATH_PRIM_MAT_FUN_CHOLESKY_DECOMPOSE_GPU_HPP
#define STAN_MATH_PRIM_MAT_FUN_CHOLESKY_DECOMPOSE_GPU_HPP

#include <stan/math/prim/mat/fun/ocl_gpu.hpp>
#include <stan/math/prim/arr/fun/matrix_gpu.hpp>
#include <stan/math/prim/mat/fun/basic_matrix_gpu.hpp>
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
      matrix_gpu A(m);
      cl::Kernel kernel_chol_block = get_kernel("cholesky_block");
      cl::Kernel kernel_left = get_kernel("cholesky_left_update");
      cl::Kernel kernel_mid = get_kernel("cholesky_mid_update");
      cl::Kernel kernel_zero = get_kernel("cholesky_zero");
      cl::CommandQueue cmd_queue = get_queue();
      try {
        cl::Context ctx = get_context();
        // Will be managed by the library core system
        int block = 64;
        int offset = 0;
        int local = 32;

        cl::Buffer buffer_V(ctx, CL_MEM_READ_WRITE,
         sizeof(T) * block * block * 4);
        cl::Buffer buffer_L(ctx, CL_MEM_READ_WRITE,
         sizeof(T) * block * A.rows() * 4);
        cl::Buffer buffer_D(ctx, CL_MEM_READ_WRITE,
         sizeof(T) * block * block * 4);

        kernel_chol_block.setArg(0, A.buffer());
        kernel_chol_block.setArg(1, offset);
        kernel_chol_block.setArg(2, A.rows());
        kernel_chol_block.setArg(3, block);
        kernel_chol_block.setArg(4, buffer_V);
        kernel_chol_block.setArg(5, buffer_D);

        kernel_left.setArg(0, buffer_L);
        kernel_left.setArg(1, A.buffer());
        kernel_left.setArg(2, buffer_V);
        kernel_left.setArg(3, offset);
        kernel_left.setArg(4, block);
        kernel_left.setArg(5, A.rows());

        kernel_mid.setArg(0, buffer_L);
        kernel_mid.setArg(1, A.buffer());
        kernel_mid.setArg(2, offset);
        kernel_mid.setArg(3, block);
        kernel_mid.setArg(4, A.rows());

        kernel_zero.setArg(0, A.buffer());
        kernel_zero.setArg(1, A.rows());

        int threadsLeft,  leftPad;
        int threadsMid,  midPad;
        while ((offset + block) < (A.rows() - block)) {
          threadsLeft = A.rows() - offset - block;
          leftPad = ((threadsLeft + local - 1) / local) * local;
          threadsMid = A.rows() - offset - block;
          midPad = ((threadsMid + local - 1) / local) * local;
          kernel_left.setArg(3, offset);
          kernel_left.setArg(6, threadsLeft);
          kernel_mid.setArg(2, offset);
          kernel_mid.setArg(5, threadsMid);
          kernel_chol_block.setArg(1, offset);
          cmd_queue.enqueueNDRangeKernel(kernel_chol_block,
           cl::NullRange, cl::NDRange(block), cl::NDRange(block));
          cmd_queue.enqueueNDRangeKernel(kernel_left,
           cl::NullRange, cl::NDRange(leftPad, leftPad),
           cl::NDRange(local,  local));
          cmd_queue.enqueueNDRangeKernel(kernel_mid,
           cl::NullRange, cl::NDRange(midPad, midPad),
           cl::NDRange(local, local));
          offset += block;
        }
        int left = A.rows() - offset;
        if (left > 0) {
          kernel_chol_block.setArg(1, offset);
          kernel_chol_block.setArg(3, left);
          cmd_queue.enqueueNDRangeKernel(kernel_chol_block,
           cl::NullRange, cl::NDRange(left), cl::NDRange(left));
        }

      cmd_queue.enqueueNDRangeKernel(kernel_zero,
       cl::NullRange, cl::NDRange(A.rows(), A.rows()), cl::NullRange);
      } catch (const cl::Error& e) {
        check_ocl_error(e);
      }
      copy_triangular_transposed(A, LOWER_TO_UPPER_TRIANGULAR);
      check_positive_definite_gpu("cholesky_decompose_gpu",
        "Matrix m", A);
      Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
       m_tmp(m.rows(), m.cols());
      copy(A, m_tmp); // NOLINT
      return m_tmp.template triangularView<Eigen::Lower>();
    }
  }
}

#endif
