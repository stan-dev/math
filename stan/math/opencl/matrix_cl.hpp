#ifndef STAN_MATH_OPENCL_MATRIX_CL_HPP
#define STAN_MATH_OPENCL_MATRIX_CL_HPP
#ifdef STAN_OPENCL
#include <stan/math/opencl/opencl_context.hpp>
#include <stan/math/opencl/kernel_cl.hpp>
#include <stan/math/opencl/constants.hpp>
#include <stan/math/opencl/cache_copy.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/scal/err/check_size_match.hpp>
#include <stan/math/prim/scal/err/domain_error.hpp>
#include <stan/math/prim/scal/meta/is_var.hpp>
#include <stan/math/opencl/kernels/copy.hpp>
#include <stan/math/opencl/kernels/sub_block.hpp>
#include <stan/math/opencl/kernels/triangular_transpose.hpp>
#include <stan/math/opencl/kernels/zeros.hpp>
#include <stan/math/rev/scal/fun/value_of_rec.hpp>
#include <CL/cl.hpp>
#include <boost/optional.hpp>
#include <iostream>
#include <string>
#include <type_traits>
#include <vector>
#include <algorithm>

/**
 *  @file stan/math/opencl/matrix_cl.hpp
 *  @brief The matrix_cl class - allocates memory space on the OpenCL device,
 *    functions for transfering matrices to and from OpenCL devices
 */
namespace stan {
namespace math {
/**
 * Represents a matrix on the OpenCL device.
 *
 * The matrix data is stored in the oclBuffer_.
 */
class matrix_cl {
 private:
  /**
   * cl::Buffer provides functionality for working with the OpenCL buffer.
   * An OpenCL buffer allocates the memory in the device that
   * is provided by the context.
   */
  cl::Buffer oclBuffer_;
  const int rows_;
  const int cols_;

 public:
  int rows() const { return rows_; }

  int cols() const { return cols_; }

  int size() const { return rows_ * cols_; }
  const cl::Buffer& buffer() const { return oclBuffer_; }
  matrix_cl() : rows_(0), cols_(0) {}

  matrix_cl(cl::Buffer& A, const int& R, const int& C) : rows_(R), cols_(C) {
    oclBuffer_ = A;
  }

  matrix_cl(const matrix_cl& A) : rows_(A.rows()), cols_(A.cols()) {
    if (A.size() == 0)
      return;
    // the context is needed to create the buffer object
    cl::Context& ctx = opencl_context.context();
    try {
      // creates a read&write object for "size" double values
      // in the provided context
      oclBuffer_ = cl::Buffer(ctx, CL_MEM_READ_WRITE, sizeof(double) * size());

      opencl_kernels::copy(cl::NDRange(rows_, cols_), A.buffer(),
                           this->buffer(), rows_, cols_);
    } catch (const cl::Error& e) {
      check_opencl_error("copy (OpenCL)->(OpenCL)", e);
    }
  }
  /**
   * Constructor for the matrix_cl that
   * only allocates the buffer on the OpenCL device.
   *
   * @param rows number of matrix rows, must be greater or equal to 0
   * @param cols number of matrix columns, must be greater or equal to 0
   *
   * @throw <code>std::system_error</code> if the
   * matrices do not have matching dimensions
   *
   */
  matrix_cl(int rows, int cols) : rows_(rows), cols_(cols) {
    if (size() > 0) {
      cl::Context& ctx = opencl_context.context();
      try {
        // creates the OpenCL buffer of the provided size
        oclBuffer_ = cl::Buffer(ctx, CL_MEM_READ_WRITE,
                                sizeof(double) * rows_ * cols_);
      } catch (const cl::Error& e) {
        check_opencl_error("matrix constructor", e);
      }
    }
  }

  /**
   * Constructor for the matrix_cl that
   * creates a copy of the Eigen matrix on the OpenCL device.
   *
   * @tparam T type of matrix
   * @tparam R rows of matrix
   * @tparam C cols of matrix
   * @param A Eigen matrix
   *
   * @throw <code>std::system_error</code> if the
   * matrices do not have matching dimensions
   */
  template <int R, int C>
  explicit matrix_cl(const Eigen::Matrix<double, R, C>& A)
      : rows_(A.rows()), cols_(A.cols()) {
    cl::Context& ctx = opencl_context.context();
    cl::CommandQueue& queue = opencl_context.queue();
    if (A.size() > 0) {
      oclBuffer_
          = cl::Buffer(ctx, CL_MEM_READ_WRITE, sizeof(double) * A.size());
      internal::cache_copy(oclBuffer_, A);
    }
  }

  matrix_cl& operator=(const matrix_cl& a) {
    check_size_match("assignment of (OpenCL) matrices", "source.rows()",
                     a.rows(), "destination.rows()", rows());
    check_size_match("assignment of (OpenCL) matrices", "source.cols()",
                     a.cols(), "destination.cols()", cols());
    oclBuffer_ = a.buffer();
    return *this;
  }

  /**
   * Constructor for the matrix_cl that
   * creates a copy of a var type Eigen matrix on the GPU.
   *
   *
   * @tparam R rows of matrix
   * @tparam C cols of matrix
   * @param A the Eigen matrix
   *
   * @throw <code>std::system_error</code> if the
   * matrices do not have matching dimensions
   */
  template <int R, int C>
  explicit matrix_cl(const Eigen::Matrix<var, R, C>& A)
      : rows_(A.rows()), cols_(A.cols()) {
    cl::Context& ctx = opencl_context.context();
    cl::CommandQueue& queue = opencl_context.queue();
    if (A.size() > 0) {
      Eigen::Matrix<double, -1, -1> L_A(value_of_rec(A));
      oclBuffer_
          = cl::Buffer(ctx, CL_MEM_READ_WRITE, sizeof(double) * L_A.size());
      queue.enqueueWriteBuffer(oclBuffer_, CL_TRUE, 0,
                               sizeof(double) * L_A.size(), L_A.data());
    }
  }

  /**
   * Stores zeros in the matrix on the OpenCL device.
   * Supports writing zeroes to the lower and upper triangular or
   * the whole matrix.
   *
   * @tparam triangular_view Specifies if zeros are assigned to
   * the entire matrix, lower triangular or upper triangular. The
   * value must be of type TriangularViewCL
   */
  template <TriangularViewCL triangular_view = TriangularViewCL::Entire>
  void zeros() {
    if (size() == 0)
      return;
    cl::CommandQueue cmdQueue = opencl_context.queue();
    try {
      opencl_kernels::zeros(cl::NDRange(this->rows(), this->cols()),
                            this->buffer(), this->rows(), this->cols(),
                            triangular_view);
    } catch (const cl::Error& e) {
      check_opencl_error("zeros", e);
    }
  }

  /**
   * Copies a lower/upper triangular of a matrix to it's upper/lower.
   *
   * @tparam triangular_map Specifies if the copy is
   * lower-to-upper or upper-to-lower triangular. The value
   * must be of type TriangularMap
   *
   * @throw <code>std::invalid_argument</code> if the matrix is not square.
   *
   */
  template <TriangularMapCL triangular_map = TriangularMapCL::LowerToUpper>
  void triangular_transpose() {
    if (size() == 0 || size() == 1) {
      return;
    }
    check_size_match("triangular_transpose ((OpenCL))",
                     "Expecting a square matrix; rows of ", "A", rows(),
                     "columns of ", "A", cols());

    cl::CommandQueue cmdQueue = opencl_context.queue();
    try {
      opencl_kernels::triangular_transpose(
          cl::NDRange(this->rows(), this->cols()), this->buffer(), this->rows(),
          this->cols(), triangular_map);
    } catch (const cl::Error& e) {
      check_opencl_error("triangular_transpose", e);
    }
  }
  /**
   * Write the context of A into
   * <code>this</code> starting at the top left of <code>this</code>
   * @param A input matrix
   * @param A_i the offset row in A
   * @param A_j the offset column in A
   * @param this_i the offset row for the matrix to be subset into
   * @param this_j the offset col for the matrix to be subset into
   * @param nrows the number of rows in the submatrix
   * @param ncols the number of columns in the submatrix
   */
  template <TriangularViewCL triangular_view = TriangularViewCL::Entire>
  void sub_block(const matrix_cl& A, int A_i, int A_j, int this_i, int this_j,
                 int nrows, int ncols) {
    if (nrows == 0 || ncols == 0) {
      return;
    }
    if ((A_i + nrows) > A.rows() || (A_j + ncols) > A.cols()
        || (this_i + nrows) > this->rows() || (this_j + ncols) > this->cols()) {
      domain_error("sub_block", "submatrix in *this", " is out of bounds", "");
    }
    cl::CommandQueue cmdQueue = opencl_context.queue();
    try {
      opencl_kernels::sub_block(cl::NDRange(nrows, ncols), A.buffer(),
                                this->buffer(), A_i, A_j, this_i, this_j, nrows,
                                ncols, A.rows(), A.cols(), this->rows(),
                                this->cols(), triangular_view);
    } catch (const cl::Error& e) {
      check_opencl_error("copy_submatrix", e);
    }
  }
};

class matrix_v_cl {
 private:
  /**
   * cl::Buffer provides functionality for working with the OpenCL buffer.
   * An OpenCL buffer allocates the memory in the device that
   * is provided by the context.
   */
  const int rows_;
  const int cols_;

 public:
  matrix_cl val_;
  matrix_cl adj_;
  int rows() const { return rows_; }

  int cols() const { return cols_; }

  int size() const { return rows_ * cols_; }

  matrix_v_cl() : rows_(0), cols_(0) {}

  /**
   * Constructor for the matrix_v_cl
   * for triangular matrices.
   *
   *
   * @param A the Eigen matrix
   * @param M The Rows and Columns of the matrix
   *
   * @throw <code>std::system_error</code> if the
   * matrices do not have matching dimensions
   */
  matrix_v_cl(vari**& A, const int& M)
      : rows_(M), cols_(M), val_(M, M), adj_(M, M) {
    cl::Context& ctx = opencl_context.context();
    cl::CommandQueue& queue = opencl_context.queue();
    if (size() > 0) {
      try {
        const int vari_size = M * M;
        // creates the OpenCL buffer to copy the Eigen
        // matrix to the OpenCL device
        std::vector<double> val_cpy;
        std::vector<double> adj_cpy;
        val_cpy.reserve(vari_size);
        adj_cpy.reserve(vari_size);
        // Make a flat version of a lower triangular
        size_t pos = 0;
        for (size_type j = 0; j < M; ++j) {
          for (size_type k = 0; k < j; ++k) {
            val_cpy.push_back(0);
            adj_cpy.push_back(0);
          }
          for (size_type i = j; i < M; ++i) {
            val_cpy.push_back(A[pos]->val_);
            adj_cpy.push_back(A[pos]->adj_);
            ++pos;
          }
        }
        printf("val size: %i \n ", val_cpy.size());
        cl::Buffer oclBuffer_val
            = cl::Buffer(ctx, CL_MEM_READ_WRITE, sizeof(double) * vari_size);
        cl::Buffer oclBuffer_adj
            = cl::Buffer(ctx, CL_MEM_READ_WRITE, sizeof(double) * vari_size);
        queue.enqueueWriteBuffer(oclBuffer_val, CL_TRUE, 0,
                                 sizeof(double) * vari_size, val_cpy.data());
        queue.enqueueWriteBuffer(oclBuffer_adj, CL_TRUE, 0,
                                 sizeof(double) * vari_size, adj_cpy.data());
        val_ = matrix_cl(oclBuffer_val, M, M);
        adj_ = matrix_cl(oclBuffer_adj, M, M);
      } catch (const cl::Error& e) {
        check_opencl_error("matrix constructor", e);
      }
    }
  }
};

}  // namespace math
}  // namespace stan

#endif
#endif
