#ifndef STAN_MATH_PRIM_ARR_MATRIX_GPU_HPP
#define STAN_MATH_PRIM_ARR_MATRIX_GPU_HPP

#include <stan/math/prim/mat/fun/ocl_gpu.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/arr/fun/value_of.hpp>
#include <stan/math/prim/mat/err/check_matching_dims.hpp>
#include <stan/math/prim/scal/err/check_size_match.hpp>
#include <iostream>
#include <string>
#include <vector>
#if defined(__APPLE__) || defined(__MACOSX)
#include <OpenCL/cl.hpp>
#else
#include <CL/cl.h>
#endif

enum triangularity {LOWER = 0, UPPER = 1, NONE = 2 };
enum copy_transposed_triangular {LOWER_TO_UPPER_TRIANGULAR = 0,
  UPPER_TO_LOWER_TRIANGULAR = 1};

/*
*  @file stan/math/prim/mat/fun/matrix_gpu.hpp
*    @brief The matrix_gpu class - allocates memory space on the GPU,
*      functions for transfering matrices to and from the GPU
*/

namespace stan {
  namespace math {

    enum triangularity {LOWER = 0, UPPER = 1, NONE = 2 };
    enum copy_transposed_triangular {LOWER_TO_UPPER_TRIANGULAR = 0,
      UPPER_TO_LOWER_TRIANGULAR = 1};

    /**
     * This class represents a matrix on the GPU. 
     * 
     * The matrix data is stored in the oclBuffer.     
     * 
     */
    class matrix_gpu {
      private:
        cl::Buffer oclBuffer;

      public:
        int rows_;
        int cols_;

        int rows() {
          return rows_;
        }

        int cols() {
          return cols_;
        }

        int size() {
          return rows_ * cols_;
        }

        cl::Buffer& buffer() {
          return oclBuffer;
        }

        // TODO(Rok): constructors with enumerator added when
        //  the matrix_gpu does not need to be READ_WRITE
        /**
         * Constructor for the matrix_gpu that 
         * only allocates the buffer on the GPU.
         * 
         * @param rows number of matrix rows
         * @param cols number of matrix columns
         * 
         * @throw <code>std::invalid_argument</code> if the 
         * matrices do not have matching dimensions
         * 
         */
        matrix_gpu(int rows,  int cols)
        : rows_(rows), cols_(cols) {
          try {
            cl::Context ctx = get_context();
            oclBuffer = cl::Buffer(ctx, CL_MEM_READ_WRITE,
             sizeof(double) * rows_ * cols_);
          } catch (const cl::Error& e) {
            check_ocl_error(e);
          }
        }

        /**
         * Constructor for the matrix_gpu that 
         * creates a copy of the Eigen matrix on the GPU.
         * 
         * 
         * @tparam T type of data in the Eigen matrix
         * @param A the Eigen matrix         
         * 
         * @throw <code>std::invalid_argument</code> if the 
         * matrices do not have matching dimensions
         * 
         */
        template<typename T, int R, int C>
        explicit matrix_gpu(const Eigen::Matrix<T, R, C> &A) {
          try {
            cl::Context ctx = get_context();
            cl::CommandQueue queue = get_queue();
            rows_ = A.rows();
            cols_ = A.cols();
            oclBuffer = cl::Buffer(ctx, CL_MEM_READ_WRITE,
             sizeof(double) * A.size());

            cl::Buffer buffer_temp = cl::Buffer(ctx, CL_MEM_READ_WRITE,
             sizeof(double) * A.size());

            queue.enqueueWriteBuffer(buffer_temp, CL_TRUE, 0,
             sizeof(double) * A.size(), A.data());

            cl::Kernel kernel = get_kernel("transpose");

            kernel.setArg(0, oclBuffer);
            kernel.setArg(1, buffer_temp);
            kernel.setArg(2, cols_);
            kernel.setArg(3, rows_);

            queue.enqueueNDRangeKernel(
             kernel,
             cl::NullRange,
             cl::NDRange(cols_, rows_),
             cl::NullRange,
             NULL,
             NULL);
          } catch (const cl::Error& e) {
            check_ocl_error(e);
          }
        }
    };

    /**
     * Copies the source Eigen matrix to 
     * the destination matrix that is stored 
     * on the GPU.
     * 
     * @tparam T type of data in the Eigen matrix
     * @param src source Eigen matrix
     * @param dst destination matrix on the GPU
     * 
     * @throw <code>std::invalid_argument</code> if the 
     * matrices do not have matching dimensions
     * 
     */
    template <typename T, int R, int C>
    void copy(const Eigen::Matrix<T, R, C> & src,
     matrix_gpu & dst) {
            check_size_match("copy (Eigen -> GPU)",
             "src.rows()", src.rows(), "dst.rows()", dst.rows());
            check_size_match("copy (Eigen -> GPU)",
             "src.cols()", src.cols(), "dst.cols()", dst.cols());

            try {
              cl::Context ctx = get_context();
              cl::CommandQueue queue = get_queue();
              cl::Buffer buffer = dst.buffer();

              cl::Buffer buffer_temp = cl::Buffer(ctx, CL_MEM_READ_WRITE,
               sizeof(T) * src.size());

              queue.enqueueWriteBuffer(buffer_temp, CL_TRUE, 0,
               sizeof(T) * dst.size(), src.data());

              cl::Kernel kernel = get_kernel("transpose");

              kernel.setArg(0, buffer);
              kernel.setArg(1, buffer_temp);
              kernel.setArg(2, src.cols());
              kernel.setArg(3, src.rows());
              queue.enqueueNDRangeKernel(
                  kernel,
                  cl::NullRange,
                  cl::NDRange(src.cols(), src.rows()),
                  cl::NullRange,
                  NULL,
                  NULL);
            } catch (const cl::Error& e) {
              check_ocl_error(e);
            }
    }

    /**
     * Copies the source matrix that is stored 
     * on the GPU to the destination Eigen 
     * matrix. 
     * 
     * @tparam T type of data in the Eigen matrix
     * @param src source matrix on the GPU
     * @param dst destination Eigen matrix
     * 
     * @throw <code>std::invalid_argument</code> if the 
     * matrices do not have matching dimensions
     * 
     */
    template <typename T, int R, int C>
    void copy(matrix_gpu & src,
     Eigen::Matrix<T, R, C> & dst) {
            check_size_match("copy (GPU -> Eigen)",
             "src.rows()", src.rows(), "dst.rows()", dst.rows());
            check_size_match("copy (GPU -> Eigen)",
             "src.cols()", src.cols(), "dst.cols()", dst.cols());
            try {
              cl::Context ctx = get_context();
              cl::CommandQueue queue = get_queue();
              cl::Buffer buffer = src.buffer();

              cl::Buffer buffer_temp = cl::Buffer(ctx, CL_MEM_READ_WRITE,
                sizeof(T) * src.size());

              cl::Kernel kernel = get_kernel("transpose");

              kernel.setArg(0, buffer_temp);
              kernel.setArg(1, buffer);
              kernel.setArg(2, src.rows());
              kernel.setArg(3, src.cols());
              queue.enqueueNDRangeKernel(
                  kernel,
                  cl::NullRange,
                  cl::NDRange(src.rows(), src.cols()),
                  cl::NullRange,
                  NULL,
                  NULL);

              queue.enqueueReadBuffer(buffer_temp,  CL_TRUE,  0,
                sizeof(T) * dst.size(),  dst.data());
            } catch (const cl::Error& e) {
              check_ocl_error(e);
            }
    }

    /**
     * Copies the source matrix to the 
     * destination matrix. Both matrices 
     * are stored on the GPU.
     * 
     * @param src source matrix
     * @param dst destination matrix
     * 
     * @throw <code>std::invalid_argument</code> if the 
     * matrices do not have matching dimensions
     * 
     */
    void copy(matrix_gpu& src,  matrix_gpu& dst) { // NOLINT
      check_size_match("copy (GPU -> GPU)",
        "src.rows()", src.rows(), "dst.rows()", dst.rows());
      check_size_match("copy (GPU -> GPU)",
        "src.cols()", src.cols(), "dst.cols()", dst.cols());
      cl::Kernel kernel = get_kernel("copy");
      cl::CommandQueue cmdQueue = get_queue();
      try {
        kernel.setArg(0, src.buffer());
        kernel.setArg(1, dst.buffer());
        kernel.setArg(2, dst.rows());
        kernel.setArg(3, dst.cols());
        cmdQueue.enqueueNDRangeKernel(
          kernel,
          cl::NullRange,
          cl::NDRange(dst.rows(), dst.cols()),
          cl::NullRange,
          NULL,
          NULL);
      } catch (const cl::Error& e) {
        check_ocl_error(e);
      }
    }

  }
}

#endif
