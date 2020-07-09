#ifndef STAN_MATH_OPENCL_REV_MATRIX_CL_HPP
#define STAN_MATH_OPENCL_REV_MATRIX_CL_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <stan/math/prim/fun/vec_concat.hpp>
#include <stan/math/rev/fun/typedefs.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <CL/cl2.hpp>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

namespace stan {
namespace math {

/** \addtogroup matrix_cl_group
 *  @{
 */

/**
 * Represents an var matrix on the OpenCL device.
 * @tparam T an arithmetic type for the type stored in the OpenCL buffer.
 */
template <typename T>
class matrix_cl<T, require_var_t<T>> {
 private:
  const int rows_{0};            // Number of rows
  const int cols_{0};            // Number of columns
  matrix_cl<double> val_{0, 0};  // holds autodiff values
  matrix_cl<double> adj_{0, 0};  // holds autodiff adjoints
  matrix_cl_view view_{matrix_cl_view::Entire};

 public:
  using Scalar = T;  // Inner type of matrix
  using type = T;    // Inner type of matrix
  /**
   * Return the number of rows
   */
  inline int rows() const { return rows_; }

  /**
   * Return the number of cols
   */
  inline int cols() const { return cols_; }

  /**
   * Return the size of the matrix
   */
  inline int size() const { return rows_ * cols_; }

  /**
   * Constant accessor value matrix
   */
  inline const matrix_cl<double>& val() const { return val_; }

  /**
   * Accessor for value matrix
   */
  inline matrix_cl<double>& val() { return val_; }

  /**
   * Constant accessor for the adjoint matrix
   */
  inline const matrix_cl<double>& adj() const { return adj_; }
  /**
   * Accessor for the adjoint matrix
   */
  inline matrix_cl<double>& adj() { return adj_; }

  matrix_cl() : rows_(0), cols_(0), val_(0, 0), adj_(0, 0) {}

  // Forward declare the methods that work in place on the matrix
  template <matrix_cl_view matrix_view = matrix_cl_view::Entire>
  inline void zeros();
  template <matrix_cl_view matrix_view = matrix_cl_view::Entire>
  inline void zeros_strict_tri();
  template <TriangularMapCL triangular_map = TriangularMapCL::LowerToUpper>
  inline void triangular_transpose();

  inline void sub_block(const matrix_cl<T, require_var_t<T>>& A, size_t A_i,
                        size_t A_j, size_t this_i, size_t this_j, size_t nrows,
                        size_t ncols);

  /**
   * Read only accessor of the view
   */
  inline const matrix_cl_view& view() const { return view_; }

  /**
   * Modify the view.
   * @param view a `matrix_cl_view` that indicates a special type of matrix.
   */
  inline void view(const matrix_cl_view& view) {
    view_ = view;
    val_.view(view);
    adj_.view(view);
  }

  /**
   * Construct a matrix from an Eigen type.
   * @tparam Mat Type derived from `Eigen::Base`
   * @param A an object derived from `Eigen::EigenBase`
   * @param partial_view `matrix_cl_view` for declaring special type.
   */
  template <typename Mat, require_eigen_st<is_var, Mat>* = nullptr>
  explicit matrix_cl(Mat&& A,
                     matrix_cl_view partial_view = matrix_cl_view::Entire)
      : rows_(A.rows()),
        cols_(A.cols()),
        val_(A.val().eval(), partial_view),
        adj_(A.adj().eval(), partial_view),
        view_(partial_view) {}

  /**
   * Construct a matrix from an `Eigen::Matrix<var, -1, -1>`
   * @param A Eigen matrix with a scalar `var` type.
   * @param partial_view `matrix_cl_view` for declaring special type.
   */
  explicit matrix_cl(const matrix_vi& A,
                     matrix_cl_view partial_view = matrix_cl_view::Entire)
      : rows_(A.rows()),
        cols_(A.cols()),
        val_(A.val().eval(), partial_view),
        adj_(A.adj().eval(), partial_view),
        view_(partial_view) {}

  /**
   * Construct a matrix from a pointer of vari pointers
   * @param A a pointer pointing to varis
   * @param R number of rows for the matrix.
   * @param C Number of columns for the matrix.
   * @param partial_view `matrix_cl_view` for declaring special type.
   */
  explicit matrix_cl(vari** A, const int& R, const int& C,
                     matrix_cl_view partial_view = matrix_cl_view::Entire)
      : rows_(R),
        cols_(C),
        view_(partial_view),
        val_(R, C, partial_view),
        adj_(R, C, partial_view) {
    cl::Context& ctx = opencl_context.context();
    cl::CommandQueue& queue = opencl_context.queue();
    cl::Event transfer_event1;
    cl::Event transfer_event2;
    if (view_ == matrix_cl_view::Entire) {
      const matrix_vi foo = Eigen::Map<const matrix_vi>(A, R, C);
      val_.buffer()
          = cl::Buffer(ctx, CL_MEM_READ_WRITE, sizeof(double) * R * C);
      adj_.buffer()
          = cl::Buffer(ctx, CL_MEM_READ_WRITE, sizeof(double) * R * C);
      queue.enqueueWriteBuffer(val_.buffer(), CL_FALSE, 0,
                               sizeof(double) * R * C, foo.val().eval().data(),
                               NULL, &transfer_event1);
      queue.enqueueWriteBuffer(adj_.buffer(), CL_FALSE, 0,
                               sizeof(double) * R * C, foo.adj().eval().data(),
                               NULL, &transfer_event2);
    } else {
      const int packed_size = R * (R + 1) / 2;
      val_.buffer()
          = cl::Buffer(ctx, CL_MEM_READ_WRITE, sizeof(double) * packed_size);
      adj_.buffer()
          = cl::Buffer(ctx, CL_MEM_READ_WRITE, sizeof(double) * packed_size);
      const matrix_vi foo = Eigen::Map<const matrix_vi>(A, 1, packed_size);
      queue.enqueueWriteBuffer(val_.buffer(), CL_FALSE, 0,
                               sizeof(double) * packed_size,
                               foo.val().eval().data(), NULL, &transfer_event1);
      queue.enqueueWriteBuffer(adj_.buffer(), CL_FALSE, 0,
                               sizeof(double) * packed_size,
                               foo.adj().eval().data(), NULL, &transfer_event2);
    }
    val_.add_write_event(transfer_event1);
    adj_.add_write_event(transfer_event2);
  }

  /**
   * Create matrix of vars from a vector.
   * @param A a standard vector holding vars
   * @param R number of rows for the matrix.
   * @param C Number of columns for the matrix.
   * @param partial_view `matrix_cl_view` for declaring special type.
   */
  explicit matrix_cl(const std::vector<var>& A, const int& R, const int& C,
                     matrix_cl_view partial_view = matrix_cl_view::Entire)
      : rows_(R),
        cols_(C),
        view_(partial_view),
        val_(R, C, partial_view),
        adj_(R, C, partial_view) {
    cl::Context& ctx = opencl_context.context();
    cl::CommandQueue& queue = opencl_context.queue();
    cl::Event transfer_event1;
    cl::Event transfer_event2;
    if (view_ == matrix_cl_view::Entire) {
      const matrix_v foo = Eigen::Map<const matrix_v>(A.data(), R, C);
      val_.buffer()
          = cl::Buffer(ctx, CL_MEM_READ_WRITE, sizeof(double) * R * C);
      adj_.buffer()
          = cl::Buffer(ctx, CL_MEM_READ_WRITE, sizeof(double) * R * C);
      queue.enqueueWriteBuffer(val_.buffer(), CL_FALSE, 0,
                               sizeof(double) * R * C, foo.val().eval().data(),
                               NULL, &transfer_event1);
      queue.enqueueWriteBuffer(adj_.buffer(), CL_FALSE, 0,
                               sizeof(double) * R * C, foo.adj().eval().data(),
                               NULL, &transfer_event2);
    } else {
      const int packed_size = R * (R + 1) / 2;
      val_.buffer()
          = cl::Buffer(ctx, CL_MEM_READ_WRITE, sizeof(double) * packed_size);
      adj_.buffer()
          = cl::Buffer(ctx, CL_MEM_READ_WRITE, sizeof(double) * packed_size);
      const matrix_v foo = Eigen::Map<const matrix_v>(A.data(), 1, packed_size);
      queue.enqueueWriteBuffer(val_.buffer(), CL_FALSE, 0,
                               sizeof(double) * packed_size,
                               foo.val().eval().data(), NULL, &transfer_event1);
      queue.enqueueWriteBuffer(adj_.buffer(), CL_FALSE, 0,
                               sizeof(double) * packed_size,
                               foo.adj().eval().data(), NULL, &transfer_event2);
    }
    val_.add_write_event(transfer_event1);
    adj_.add_write_event(transfer_event2);
  }

  /**
   * Initialize a var matrix with size rows and columns
   */
  explicit matrix_cl(const int& rows, const int& cols,
                     matrix_cl_view partial_view = matrix_cl_view::Entire)
      : rows_(rows),
        cols_(cols),
        val_(rows, cols),
        adj_(rows, cols),
        view_(partial_view) {}

  auto& operator=(const matrix_cl<var>& A) {
    val_ = A.val();
    adj_ = A.adj();
    return *this;
  }
};

/** @}*/

}  // namespace math
}  // namespace stan

#endif
#endif
