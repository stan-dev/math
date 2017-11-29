#ifndef STAN_MATH_PRIM_ARR_MATRIX_GPU_HPP
#define STAN_MATH_PRIM_ARR_MATRIX_GPU_HPP

#include <stan/math/prim/mat/fun/ocl_gpu.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/arr/fun/value_of.hpp>
#include <stan/math/prim/mat/err/check_matching_dims.hpp>
#include <stan/math/prim/scal/err/check_size_match.hpp>
#include <CL/cl.hpp>
#include <iostream>
#include <string>
#include <vector>

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
     * The matrix data is stored in the oclBuffer_.     
     * 
     */
    class matrix_gpu {
      private:
        cl::Buffer oclBuffer_;
        const int rows_;
        const int cols_;

      public:
        int rows() const {
          return rows_;
        }

        int cols() const  {
          return cols_;
        }

        int size() const {
          return rows_ * cols_;
        }

        const cl::Buffer& buffer() const {
          return oclBuffer_;
        }

        matrix_gpu()
         : rows_(0), cols_(0) {
        }

        matrix_gpu(matrix_gpu& a)
         : rows_(a.rows()), cols_(a.cols()) {
        }

        matrix_gpu(const matrix_gpu& a)
         : rows_(a.rows()), cols_(a.cols()) {
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
          if ( size() > 0 ) {
            cl::Context& ctx = get_context();
            try {
              oclBuffer_ = cl::Buffer(ctx, CL_MEM_READ_WRITE,
                 sizeof(double) * rows_ * cols_);
            } catch (const cl::Error& e) {
              check_ocl_error("matrix constructor", e);
            }
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
        explicit matrix_gpu(const Eigen::Matrix<T, R, C> &A)
        : rows_(A.rows()), cols_(A.cols()) {
          if ( size() > 0 ) {
            cl::Context& ctx = get_context();
            cl::CommandQueue& queue = get_queue();
            try {
              oclBuffer_ = cl::Buffer(ctx, CL_MEM_READ_WRITE,
               sizeof(double) * A.size());

              queue.enqueueWriteBuffer(oclBuffer_, CL_TRUE, 0,
               sizeof(T) * A.size(), A.data());
            } catch (const cl::Error& e) {
              check_ocl_error("matrix constructor", e);
            }
          }
        }

        matrix_gpu& operator= (const matrix_gpu& a){
          check_size_match("assignment of GPU matrices",
           "source.rows()", a.rows(), "destination.rows()", rows());
          check_size_match("assignment of GPU matrices",
           "source.cols()", a.cols(), "destination.cols()", cols());
          oclBuffer_ = a.buffer();
          return *this;
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
    void copy(const Eigen::Matrix<T, R, C>& src,
     matrix_gpu& dst) {
      check_size_match("copy (Eigen -> GPU)",
       "src.rows()", src.rows(), "dst.rows()", dst.rows());
      check_size_match("copy (Eigen -> GPU)",
       "src.cols()", src.cols(), "dst.cols()", dst.cols());
      if ( src.size() > 0 ) {
        cl::CommandQueue queue = get_queue();
        try {
          queue.enqueueWriteBuffer(dst.buffer(), CL_TRUE, 0,
           sizeof(T) * dst.size(), src.data());
        } catch (const cl::Error& e) {
          check_ocl_error("copy Eigen->GPU", e);
        }
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
      if ( src.size() > 0 ) {
        cl::CommandQueue queue = get_queue();
        try {
          queue.enqueueReadBuffer(src.buffer(),  CL_TRUE,  0,
            sizeof(T) * dst.size(),  dst.data());
        } catch (const cl::Error& e) {
          check_ocl_error("copy GPU->Eigen", e);
        }
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
      if ( src.size() > 0 ) {
        cl::Kernel kernel = get_kernel("copy");
        cl::CommandQueue& cmdQueue = get_queue();
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
          check_ocl_error("copy GPU->GPU", e);
        }
      }
    }

  }
}

#endif
