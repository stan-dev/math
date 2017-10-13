#ifndef STAN_MATH_PRIM_ARR_MATRIX_GPU_HPP
#define STAN_MATH_PRIM_ARR_MATRIX_GPU_HPP

#include <stan/math/prim/mat/fun/ocl_gpu.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/arr/fun/value_of.hpp>
#include <iostream>
#include <string>
#if defined(__APPLE__) || defined(__MACOSX)
#include <OpenCL/cl.hpp>
#else
#include <CL/cl.hpp>
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

    class matrix_gpu {
      private:

        cl::Buffer oclBuffer;

      public:

        int rows;
        int cols;
        cl::Buffer buffer() {
          return oclBuffer;
        }

        //TODO: constructors with enumerator added when
        // the matrix_gpu does not need to be READ_WRITE
        matrix_gpu(int rows,  int cols) {
          try {
            cl::Context ctx = stan::math::get_context();
            oclBuffer = cl::Buffer(ctx, CL_MEM_READ_WRITE,
             sizeof(double) * rows * cols);
            this->rows = rows;
            this->cols = cols;
          } catch (cl::Error& e) {
            check_ocl_error(e);
          }
        };

        template <typename T>
        matrix_gpu(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &A) {
          try {
            cl::Context ctx = stan::math::get_context();
            cl::CommandQueue queue = stan::math::get_queue();
            this->rows = A.rows();
            this->cols = A.cols();
            oclBuffer = cl::Buffer(ctx, CL_MEM_READ_WRITE, sizeof(double) *
             A.rows() * A.cols());

            double*  Atemp =  new double[rows * cols];
            for(int j = 0; j < cols; j++) {
              for(int i = 0; i < rows; i++) {
      //            std::cout << "value_of(A(i,j)): " << value_of(A(i,j)) << "\n";
                  Atemp[i * cols + j] = value_of(A(i, j));
      //             std::cout << "Atemp[i * cols + j]: " << Atemp[i * cols + j] << "\n";
              }
            }
            queue.enqueueWriteBuffer(oclBuffer, CL_TRUE, 0,
             sizeof(double) * A.rows() * A.cols(), Atemp);
            delete[] Atemp;
          } catch (cl::Error& e) {
            check_ocl_error(e);
          }
        };
    };

    template <typename T>
    void copy(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> & src,
     stan::math::matrix_gpu & dst) {
            if (src.rows() != dst.rows || src.cols() != dst.cols)
            {
              std::cout << "Copy (Eigen -> GPU) :Eigen and matrix_gpu sizes do no match!" <<
                std::endl;
            }
            int rows = src.rows();
            int cols = src.cols();
            try {
              cl::Context ctx = stan::math::get_context();
              cl::CommandQueue queue = stan::math::get_queue();
              cl::Buffer buffer = dst.buffer();
              double*  Atemp =  new double[rows * cols];
              for(int j = 0; j < cols; j++) {
                for(int i = 0; i < rows; i++) {
                  Atemp[i * cols + j] = value_of(src(i, j));
                }
              }
              queue.enqueueWriteBuffer(buffer, CL_TRUE, 0,
               sizeof(double) * rows * cols, Atemp);
              delete[] Atemp;

            } catch (cl::Error& e) {
              check_ocl_error(e);
            }
    }

    template <typename T>
    void copy(stan::math::matrix_gpu & src,
     Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> & dst) {
            if (dst.rows() != src.rows) {
              std::cout << "Copy (GPU -> Eigen): Eigen and matrix_gpu row lengths do no match!" <<
               std::endl;
            }
            if (dst.cols() != src.cols) {
              std::cout << "Copy (GPU -> Eigen): Eigen and stanmathCL matrix_gpu col lengths do no match!" <<
               std::endl;
            }
            int rows = dst.rows();
            int cols = dst.cols();
            try {
              double*  Btemp =  new double[rows * cols];
              cl::Context ctx = stan::math::get_context();
              cl::CommandQueue queue = stan::math::get_queue();
              cl::Buffer buffer = src.buffer();
              queue.enqueueReadBuffer(buffer,  CL_TRUE,  0,
                sizeof(double) * dst.size(),  Btemp);
              for(int j = 0; j < cols; j++) {
                for(int i = 0; i < rows; i++) {
                  dst(i, j) = Btemp[i * cols + j];
                }
              }
              delete[] Btemp;
            } catch (cl::Error& e) {
              check_ocl_error(e);
            }
    }

    void copy(stan::math::matrix_gpu & src,  stan::math::matrix_gpu & dst) {
      if (src.rows != dst.rows || src.cols != dst.cols) {
        app_error("The dimensions of the input and output matrices should match.");
      }
      cl::Kernel kernel = stan::math::get_kernel("copy");
      cl::CommandQueue cmdQueue = stan::math::get_queue();
      try {
        kernel.setArg(0, src.buffer());
        kernel.setArg(1, dst.buffer());
        kernel.setArg(2, dst.rows);
        kernel.setArg(3, dst.cols);
        cmdQueue.enqueueNDRangeKernel(
          kernel,
          cl::NullRange,
          cl::NDRange(dst.rows, dst.cols),
          cl::NullRange,
          NULL,
          NULL);
      } catch (cl::Error& e) {
        check_ocl_error(e);
      }
    }

  }
}

#endif
