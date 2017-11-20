#ifndef STAN_MATH_PRIM_MAT_ERR_CHECK_GPU_HPP
#define STAN_MATH_PRIM_MAT_ERR_CHECK_GPU_HPP

#include <stan/math/prim/mat/fun/ocl_gpu.hpp>
#include <stan/math/prim/arr/fun/matrix_gpu.hpp>
#include <stan/math/prim/scal/err/domain_error.hpp>
#include <stan/math/prim/scal/err/check_size_match.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <sstream>
#include <string>

namespace stan {
  namespace math {
     /**
     * Check if the specified matrix on the GPU is square.
     *
     * This check allows 0x0 matrices.
     *
     *
     * @param function Function name (for error messages)
     * @param name Variable name (for error messages)
     * @param y Matrix on the GPU to test
     *
     * @throw <code>std::invalid_argument</code> if the matrix
     *    is not square
     */
    inline void
    check_square(const std::string& function,
                 const std::string& name,
                 matrix_gpu& y) {
      check_size_match(function,
                       "Expecting a square matrix; rows of ", name, y.rows(),
                       "columns of ", name, y.cols());
    }
    /**
     * Check if the specified matrix on the GPU has NaN values
     *
     * @param function Function name (for error messages)
     * @param name Variable name (for error messages)
     * @param y Matrix to test
     *
     * @throw <code>std::domain_error</code> if
     *    any element of the matrix is <code>NaN</code>.
     */
    inline void
    check_positive_definite_gpu(const std::string& function,
                           const std::string& name,
                  matrix_gpu& y) {
      if (y.size() == 0) return;

      cl::Kernel kernel_check_nan = get_kernel("check_nan");
      cl::Kernel kernel_check_diagonal_zeros =
       get_kernel("check_diagonal_zeros");
      cl::CommandQueue cmd_queue = get_queue();

      try {
        cl::Context& ctx = get_context();
        int nan_flag = 0;
        int diag_zeros_flag = 0;
        cl::Buffer buffer_flag(ctx, CL_MEM_READ_WRITE,
         sizeof(int));

        cmd_queue.enqueueWriteBuffer(buffer_flag, CL_TRUE, 0,
         sizeof(int), &nan_flag);

        kernel_check_nan.setArg(0, y.buffer());
        kernel_check_nan.setArg(1, y.rows());
        kernel_check_nan.setArg(2, y.cols());
        kernel_check_nan.setArg(3, buffer_flag);

        cmd_queue.enqueueNDRangeKernel(kernel_check_nan,
         cl::NullRange, cl::NDRange(y.rows(), y.cols()),  cl::NullRange);

        cmd_queue.enqueueReadBuffer(buffer_flag, CL_TRUE, 0,
         sizeof(int), &nan_flag);

        kernel_check_diagonal_zeros.setArg(0, y.buffer());
        kernel_check_diagonal_zeros.setArg(1, y.rows());
        kernel_check_diagonal_zeros.setArg(2, y.cols());
        kernel_check_diagonal_zeros.setArg(3, buffer_flag);

        cmd_queue.enqueueNDRangeKernel(kernel_check_diagonal_zeros,
         cl::NullRange, cl::NDRange(y.rows(), y.cols()),  cl::NullRange);

        cmd_queue.enqueueReadBuffer(buffer_flag, CL_TRUE, 0,
         sizeof(int), &diag_zeros_flag);

        //  if NaN values were found in the matrix
        if (nan_flag || diag_zeros_flag) {
          domain_error(function, name, "is not positive definite", "");
        }
      } catch (const cl::Error& e) {
        check_ocl_error("pos_def_check", e);
      }
    }
    /**
     * Check if the specified matrix on the GPU has NaN values
     *
     * @param function Function name (for error messages)
     * @param name Variable name (for error messages)
     * @param y Matrix to test
     *
     * @throw <code>std::domain_error</code> if
     *    any element of the matrix is <code>NaN</code>.
     */
    inline void
    check_nan_gpu(const std::string& function,
                           const std::string& name,
                  matrix_gpu& y) {
      if (y.size() == 0) return;

      cl::Kernel kernel_check_nan = get_kernel("check_nan");
      cl::CommandQueue cmd_queue = get_queue();

      try {
        cl::Context& ctx = get_context();
        int nan_flag = 0;
        cl::Buffer buffer_nan_flag(ctx, CL_MEM_READ_WRITE,
         sizeof(int));

        cmd_queue.enqueueWriteBuffer(buffer_nan_flag, CL_TRUE, 0,
         sizeof(int), &nan_flag);

        kernel_check_nan.setArg(0, y.buffer());
        kernel_check_nan.setArg(1, y.rows());
        kernel_check_nan.setArg(2, y.cols());
        kernel_check_nan.setArg(3, buffer_nan_flag);

        cmd_queue.enqueueNDRangeKernel(kernel_check_nan,
         cl::NullRange, cl::NDRange(y.rows(), y.cols()),  cl::NullRange);

        cmd_queue.enqueueReadBuffer(buffer_nan_flag, CL_TRUE, 0,
         sizeof(int), &nan_flag);
        //  if NaN values were found in the matrix
        if (nan_flag) {
          domain_error(function, name, "has NaN values", "");
        }
      } catch (const cl::Error& e) {
        check_ocl_error("nan_check", e);
      }
    }
    /**
     * Check if the specified matrix on the GPU is symmetric
     *
     * @param function Function name (for error messages)
     * @param name Variable name (for error messages)
     * @param y Matrix to test
     *
     * @throw <code>std::domain_error</code> if
     *    the matrix is not symmetric.
     */
    inline void
    check_symmetric_gpu(const std::string& function,
                           const std::string& name,
                  matrix_gpu& y) {
      if (y.size() == 0) return;
      check_square(function, name, y);
      cl::Kernel kernel_check_symmetric = get_kernel("check_symmetric");
      cl::CommandQueue cmd_queue = get_queue();

      try {
        cl::Context& ctx = get_context();
        int symmetric_flag = 0;
        cl::Buffer buffer_symmetric_flag(ctx, CL_MEM_READ_WRITE,
         sizeof(int));

        cmd_queue.enqueueWriteBuffer(buffer_symmetric_flag, CL_TRUE, 0,
         sizeof(int), &symmetric_flag);

        kernel_check_symmetric.setArg(0, y.buffer());
        kernel_check_symmetric.setArg(1, y.rows());
        kernel_check_symmetric.setArg(2, y.cols());
        kernel_check_symmetric.setArg(3, buffer_symmetric_flag);
        kernel_check_symmetric.setArg(4, CONSTRAINT_TOLERANCE);

        cmd_queue.enqueueNDRangeKernel(kernel_check_symmetric,
         cl::NullRange, cl::NDRange(y.rows(), y.cols()),  cl::NullRange);

        cmd_queue.enqueueReadBuffer(buffer_symmetric_flag, CL_TRUE, 0,
         sizeof(int), &symmetric_flag);
        //  if the matrix is not symmetric
        if (symmetric_flag) {
          domain_error(function, name, "is not symmetric", "");
        }
      } catch (const cl::Error& e) {
        check_ocl_error("symmetric_check", e);
      }
    }
    /**
     * Check if the specified matrix on the GPU has zeros on the diagonal
     *
     * @param function Function name (for error messages)
     * @param name Variable name (for error messages)
     * @param y Matrix to test
     *
     * @throw <code>std::domain_error</code> if
     *    any diagonal element of the matrix is zero.
     */
    inline void
    check_diagonal_zeros(const std::string& function,
                           const std::string& name,
                  matrix_gpu& y) {
      if (y.size() == 0) return;

      cl::Kernel kernel_check_diagonal_zeros =
       get_kernel("check_diagonal_zeros");
      cl::CommandQueue cmd_queue = get_queue();

      try {
        cl::Context ctx = get_context();

        int flag = 0;

        cl::Buffer buffer_flag(ctx, CL_MEM_READ_WRITE,
         sizeof(int));

        cmd_queue.enqueueWriteBuffer(buffer_flag, CL_TRUE, 0,
         sizeof(int), &flag);

        kernel_check_diagonal_zeros.setArg(0, y.buffer());
        kernel_check_diagonal_zeros.setArg(1, y.rows());
        kernel_check_diagonal_zeros.setArg(2, y.cols());
        kernel_check_diagonal_zeros.setArg(3, buffer_flag);

        cmd_queue.enqueueNDRangeKernel(kernel_check_diagonal_zeros,
         cl::NullRange, cl::NDRange(y.rows(), y.cols()),  cl::NullRange);

        cmd_queue.enqueueReadBuffer(buffer_flag, CL_TRUE, 0,
         sizeof(int), &flag);

        //  if zeros were found on the diagonal
        if (flag) {
          domain_error(function, name, "has zeros on the diagonal.", "");
        }
      } catch (const cl::Error& e) {
        check_ocl_error("diag_zeros_check", e);
      }
    }



    /**
     * Check if the two matrices on the GPU are of the same size.
     *
     * This function checks the runtime sizes only.
     *     
     *
     * @param function Function name (for error messages)
     * @param name1 Variable name for the first matrix (for error messages)
     * @param y1 First matrix
     * @param name2 Variable name for the second matrix (for error messages)
     * @param y2 Second matrix
     *
     * @throw <code>std::invalid_argument</code>
     * if the dimensions of the matrices do not match
     */
    inline void check_matching_dims(const std::string& function,
                                    const std::string& name1,
                                    matrix_gpu& y1,
                                    const std::string& name2,
                                    matrix_gpu& y2) {
      check_size_match(function,
                       "Rows of ", name1, y1.rows(),
                       "rows of ", name2, y2.rows());
      check_size_match(function,
                       "Columns of ", name1, y1.cols(),
                       "columns of ", name2, y2.cols());
    }

  }
}
#endif
