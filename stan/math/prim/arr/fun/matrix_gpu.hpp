#ifndef STAN_MATH_PRIM_ARR_MATRIX_GPU_HPP
#define STAN_MATH_PRIM_ARR_MATRIX_GPU_HPP

#include <stan/math/prim/mat/fun/ocl.hpp>
#include <iostream>
#include <string>
#if defined(__APPLE__) || defined(__MACOSX)
#include <OpenCL/cl.hpp>
#else
#include <CL/cl.hpp>
#endif
#include <stan/math/prim/mat/fun/Eigen.hpp>

enum triangularity { LOWER = 0,  UPPER = 1,  NONE = 2 };
enum copy_transposed_triangular { LOWER_TO_UPPER_TRIANGULAR = 0,
  UPPER_TO_LOWER_TRIANGULAR = 1 };

/*
*  @file stan/math/prim/mat/fun/matrix_gpu.hpp
*    @brief The matrix_gpu class - allocates memory space on the GPU,
*      functions for transfering matrices to and from the GPU
*/

namespace stan {
  namespace math {

    class matrix_gpu{
      private:

        cl::Buffer oclBuffer;

      public:

        int M;
        int N;
        cl::Buffer buffer() {
          return oclBuffer;
        }

        //TODO: constructors with enumerator added when
        // the matrix_gpu does not need to be READ_WRITE
        matrix_gpu(int m,  int n) {
          try {
            cl::Context ctx = stan::math::get_context();
            oclBuffer = cl::Buffer(ctx, CL_MEM_READ_WRITE,
             sizeof(T) * m * n);
            this->M = m;
            this->N = n;
          } catch (cl::Error& e) {
            check_ocl_error(e);
          }
        };

        template <typename T>
        matrix_gpu(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> A) {
          try {
            cl::Context ctx = stan::math::get_context();
            cl::CommandQueue queue = stan::math::get_queue();
            oclBuffer = cl::Buffer(ctx, CL_MEM_READ_WRITE, sizeof(T) *
             A.rows() * A.cols());
            this->M = A.rows();
            this->N = A.cols();
            T*  Atemp =  new T[M * N];
            for(int i = 0; i < M; i++) {
              for(int j = 0; j < N; j++) {
                  Atemp[i * N + j] = A(i, j);
              }
            }
            queue.enqueueWriteBuffer(oclBuffer, CL_TRUE, 0,
             sizeof(T) * M * N, Atemp);
            delete Atemp;
          } catch (cl::Error& e) {
            check_ocl_error(e);
          }
        };
    };

    template <typename T>
    void copy(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> & src,
     stan::math::matrix_gpu & dst) {
            if(src.rows()! = dst.M || src.cols()! = dst.N)
            {
              std::cout << "Eigen and stanmathCL matrix_gpu sizes do no match!" <<
                std::endl;
            }
            int M = src.rows();
            int N = src.cols();
            try {
              cl::Context ctx = stan::math::get_context();
              cl::CommandQueue queue = stan::math::get_queue();
              cl::Buffer buffer = dst.buffer();
              T*  Atemp =  new T[M * N];
              for(int i = 0; i < M; i++) {
                for(int j = 0; j < N; j++) {
                  Atemp[i * N + j] = src(i, j);
                }
              }
              queue.enqueueWriteBuffer(buffer, CL_TRUE, 0,
               sizeof(T) * M * N, Atemp);
              delete Atemp;

            } catch (cl::Error& e) {
              check_ocl_error(e);
            }
    }

    template <typename T>
    void copy(stan::math::matrix_gpu & src,
     Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> & dst) {
            if(dst.rows()! = src.M) {
              std::cout << "Eigen and stanmathCL matrix_gpu sizes do no match!" <<
               std::endl;
            }
            if(dst.cols()! = src.N) {
              std::cout << "Eigen and stanmathCL matrix_gpu sizes do no match!" <<
               std::endl;
            }
            int M = dst.rows();
            int N = dst.cols();
            try {
              T*  Btemp =  new T[M * N];
              cl::Context ctx = stan::math::get_context();
              cl::CommandQueue queue = stan::math::get_queue();
              cl::Buffer buffer = src.buffer();
              queue.enqueueReadBuffer(buffer,  CL_TRUE,  0,
                sizeof(T) * M * N,  Btemp);
              for(int i = 0; i < M; i++) {
                for(int j = 0; j < N; j++) {
                  dst(i, j) = Btemp[i * N + j];
                }
              }
              delete Btemp;
            } catch (cl::Error& e) {
              check_ocl_error(e);
            }
    }

    void copy(stan::math::matrix_gpu & src,  stan::math::matrix_gpu & dst) {
      if(!(src.M =  = dst.M && src.N =  = dst.N)) {
        app_error("output matrix_gpu dimensions in matrix_gpu copy!\n"
         "\nThe dimensions of the input and output matrices should match.");
      }
      cl::Kernel kernel = stan::math::get_kernel("copy");
      cl::CommandQueue cmdQueue = stan::math::get_queue();
      try {
        kernel.setArg(0, src.buffer());
        kernel.setArg(1, dst.buffer());
        kernel.setArg(2, dst.M);
        kernel.setArg(3, dst.N);
        cmdQueue.enqueueNDRangeKernel(
          kernel,
          cl::NullRange,
          cl::NDRange(dst.M, dst.N),
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
