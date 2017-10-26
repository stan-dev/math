#ifndef STAN_MATH_PRIM_MAT_ERR_CHECK_GPU_HPP
#define STAN_MATH_PRIM_MAT_ERR_CHECK_GPU_HPP

#include <stan/math/prim/mat/fun/ocl_gpu.hpp>
#include <stan/math/prim/arr/fun/matrix_gpu.hpp>
#include <stan/math/prim/scal/err/domain_error.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <sstream>
#include <string>

namespace stan {
  namespace math {

    /**
     * Check if the specified matrix_gpu has NaN values
     *
     * @param A matrix to check for NaN values
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
        cl::Context ctx = get_context();
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
        check_ocl_error(e);
      }
    }

    /**
     * Check if the specified matrix_gpu has zeros on the diagonal
     *
     * @param A matrix to check for zeros on the diagonal
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
        check_ocl_error(e);
      }
    }
  }
}
#endif
