#ifndef STAN_MATH_PRIM_ARR_MATRIX_GPU_HPP
#define STAN_MATH_PRIM_ARR_MATRIX_GPU_HPP

#include <stan/math/prim/mat/fun/ocl_gpu.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/arr/fun/value_of.hpp>
#include <iostream>
#include <string>
#include <vector>
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

    enum triangularity {LOWER = 0, UPPER = 1, NONE = 2 };
    enum copy_transposed_triangular {LOWER_TO_UPPER_TRIANGULAR = 0,
      UPPER_TO_LOWER_TRIANGULAR = 1};

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

        template <typename T>
        explicit matrix_gpu(const Eigen::Matrix<T,
         Eigen::Dynamic, Eigen::Dynamic> &A) {
          try {
            cl::Context ctx = get_context();
            cl::CommandQueue queue = get_queue();
            rows_ = A.rows();
            cols_ = A.cols();
            oclBuffer = cl::Buffer(ctx, CL_MEM_READ_WRITE,
             sizeof(double) * A.size());

            std::vector<double> Atemp(A.size());
            for (int j = 0; j < cols_; j++) {
              for (int i = 0; i < rows_; i++) {
                  Atemp[i * rows_ + j] = value_of(A(i, j));
              }
            }
            queue.enqueueWriteBuffer(oclBuffer, CL_TRUE, 0,
             sizeof(double) * A.size(), &Atemp[0]);
          } catch (const cl::Error& e) {
            check_ocl_error(e);
          }
        }
    };

    void copy(const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> & src,
     matrix_gpu & dst) {
            if (src.rows() != dst.rows() || src.cols() != dst.cols()) {
              std::cout << "Copy (Eigen -> GPU) : CPU and GPU matrix "
               "sizes do no match!" <<
                std::endl;
            }

            try {
              cl::Context ctx = get_context();
              cl::CommandQueue queue = get_queue();
              cl::Buffer buffer = dst.buffer();

              std::vector<double> Atemp(src.size());
              for (int j = 0; j < src.cols(); j++) {
                for (int i = 0; i < src.rows(); i++) {
                  Atemp[i * src.rows() + j] = value_of(src(i, j));
                }
              }


              queue.enqueueWriteBuffer(buffer, CL_TRUE, 0,
               sizeof(double) * dst.size(), &Atemp[0]);
            } catch (const cl::Error& e) {
              check_ocl_error(e);
            }
    }

    template <typename T>
    void copy(matrix_gpu & src,
     Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> & dst) {
            if (dst.rows() != src.rows()) {
              std::cout << "Copy (GPU -> Eigen): GPU and CPU matrics "
               "row lengths do no match!" <<
               std::endl;
            }
            if (dst.cols() != src.cols()) {
              std::cout << "Copy (GPU -> Eigen): GPU and CPU matrices "
               "col lengths do not match!" <<
               std::endl;
            }
            int rows = dst.rows();
            int cols = dst.cols();
            try {
              std::vector<double> Btemp(dst.size());
              cl::Context ctx = get_context();
              cl::CommandQueue queue = get_queue();
              cl::Buffer buffer = src.buffer();
              queue.enqueueReadBuffer(buffer,  CL_TRUE,  0,
                sizeof(double) * dst.size(),  &Btemp[0]);
              for (int j = 0; j < cols; j++) {
                for (int i = 0; i < rows; i++) {
                  dst(i, j) = Btemp[i * rows + j];
                }
              }
            } catch (const cl::Error& e) {
              check_ocl_error(e);
            }
    }

    void copy(matrix_gpu & src,
     Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> & dst) {
            if (dst.rows() != src.rows()) {
              std::cout << "Copy (GPU -> Eigen): Eigen and matrix_gpu"
              " row lengths do no match!" <<
               std::endl;
            }
            if (dst.cols() != src.cols()) {
              std::cout << "Copy (GPU -> Eigen): Eigen and stanmathCL "
               "matrix_gpu col lengths do no match!" <<
               std::endl;
            }

            try {
              std::vector<double> Btemp(dst.size());
              cl::Context ctx = get_context();
              cl::CommandQueue queue = get_queue();
              cl::Buffer buffer = src.buffer();
              queue.enqueueReadBuffer(buffer,  CL_TRUE,  0,
                sizeof(double) * dst.size(),  &Btemp[0]);
              for (int j = 0; j < dst.cols(); j++) {
                for (int i = 0; i < dst.rows(); i++) {
                  dst.coeffRef(i, j) = Btemp[i * dst.rows() + j];
                }
              }
            } catch (const cl::Error& e) {
              check_ocl_error(e);
            }
    }

    void copy(matrix_gpu & src,  matrix_gpu & dst) { // NOLINT
      if (src.rows() != dst.rows() || src.cols() != dst.cols()) {
        app_error("The dimensions of the input and output "
        "matrices should match.");
      }
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
