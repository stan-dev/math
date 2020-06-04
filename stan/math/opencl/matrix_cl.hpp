#ifndef STAN_MATH_OPENCL_MATRIX_CL_HPP
#define STAN_MATH_OPENCL_MATRIX_CL_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/vec_concat.hpp>
#include <stan/math/opencl/opencl_context.hpp>
#include <stan/math/opencl/matrix_cl_view.hpp>
#include <stan/math/opencl/is_matrix_cl.hpp>
#include <stan/math/opencl/err/check_opencl.hpp>
#include <stan/math/opencl/kernel_generator/is_kernel_expression.hpp>
#include <CL/cl2.hpp>
#include <algorithm>
#include <iostream>
#include <string>
#include <type_traits>
#include <vector>

/** \ingroup opencl
 *  \defgroup matrix_cl_group Matrix
 * The matrix_cl class - allocates memory space on the OpenCL device. Operations
 * on `matrix_cl` types are executed lazily via the kernel generator
 * and async routines.
 */
namespace stan {
namespace math {

/** \addtogroup matrix_cl_group
 *  @{
 */

/**
 * Represents an arithmetic matrix on the OpenCL device.
 * @tparam T an arithmetic type for the type stored in the OpenCL buffer.
 */
template <typename T>
class matrix_cl<T, require_arithmetic_t<T>> {
 private:
  cl::Buffer buffer_cl_;  // Holds the allocated memory on the device
  int rows_{0};           // Number of rows.
  int cols_{0};           // Number of columns.
  // Holds info on if matrix is a special type
  matrix_cl_view view_{matrix_cl_view::Entire};
  mutable std::vector<cl::Event> write_events_;  // Tracks write jobs
  mutable std::vector<cl::Event> read_events_;   // Tracks reads

 public:
  using Scalar = T;  // Underlying type of the matrix
  using type = T;    // Underlying type of the matrix
  // Forward declare the methods that work in place on the matrix
  template <matrix_cl_view matrix_view = matrix_cl_view::Entire>
  inline void zeros();
  template <matrix_cl_view matrix_view = matrix_cl_view::Entire>
  inline void zeros_strict_tri();
  template <TriangularMapCL triangular_map = TriangularMapCL::LowerToUpper>
  inline void triangular_transpose();

  inline void sub_block(const matrix_cl<T, require_arithmetic_t<T>>& A,
                        size_t A_i, size_t A_j, size_t this_i, size_t this_j,
                        size_t nrows, size_t ncols);
  int rows() const { return rows_; }

  int cols() const { return cols_; }

  int size() const { return rows_ * cols_; }

  const matrix_cl_view& view() const { return view_; }

  void view(const matrix_cl_view& view) { view_ = view; }

  /**
   * Clear the write events from the event stacks.
   */
  inline void clear_write_events() const {
    write_events_.clear();
    return;
  }

  /**
   * Clear the read events from the event stacks.
   */
  inline void clear_read_events() const {
    read_events_.clear();
    return;
  }

  /**
   * Clear the write events from the event stacks.
   */
  inline void clear_read_write_events() const {
    read_events_.clear();
    write_events_.clear();
    return;
  }

  /**
   * Get the events from the event stacks.
   * @return The write event stack.
   */
  inline const std::vector<cl::Event>& write_events() const {
    return write_events_;
  }

  /**
   * Get the events from the event stacks.
   * @return The read/write event stack.
   */
  inline const std::vector<cl::Event>& read_events() const {
    return read_events_;
  }

  /**
   * Get the events from the event stacks.
   * @return The read/write event stack.
   */
  inline const std::vector<cl::Event> read_write_events() const {
    return vec_concat(this->read_events(), this->write_events());
  }

  /**
   * Add an event to the read event stack.
   * @param new_event The event to be pushed on the event stack.
   */
  inline void add_read_event(cl::Event new_event) const {
    this->read_events_.push_back(new_event);
  }

  /**
   * Add an event to the write event stack.
   * @param new_event The event to be pushed on the event stack.
   */
  inline void add_write_event(cl::Event new_event) const {
    this->write_events_.push_back(new_event);
  }

  /**
   * Add an event to the read/write event stack.
   * @param new_event The event to be pushed on the event stack.
   */
  inline void add_read_write_event(cl::Event new_event) const {
    this->read_events_.push_back(new_event);
    this->write_events_.push_back(new_event);
  }

  /**
   * Waits for the write events and clears the read event stack.
   */
  inline void wait_for_write_events() const {
    for (cl::Event e : write_events_) {
      e.wait();
    }
    write_events_.clear();
  }

  /**
   * Waits for the read events and clears the read event stack.
   */
  inline void wait_for_read_events() const {
    for (cl::Event e : read_events_) {
      e.wait();
    }
    read_events_.clear();
  }

  /**
   * Waits for read and write events to finish and clears the read, write, and
   * read/write event stacks.
   */
  inline void wait_for_read_write_events() const {
    wait_for_read_events();
    wait_for_write_events();
  }

  const cl::Buffer& buffer() const { return buffer_cl_; }
  cl::Buffer& buffer() { return buffer_cl_; }

  matrix_cl() {}
  /**
   * Construct a matrix_cl<T> from an existing cl::Buffer object. The matrix
   * directly uses given buffer - no copying is done.
   *
   * @param A the cl::Buffer object to construct the matrix from
   * @param R number of rows
   * @param C number of columns
   * @param partial_view view of the matrix
   */
  matrix_cl(cl::Buffer& A, const int R, const int C,
            matrix_cl_view partial_view = matrix_cl_view::Entire)
      : buffer_cl_(A), rows_(R), cols_(C), view_(partial_view) {}

  /**
   * Copy constructor.
   * @param A matrix_cl to copy
   */
  matrix_cl(const matrix_cl<T>& A)
      : rows_(A.rows()), cols_(A.cols()), view_(A.view()) {
    if (A.size() == 0) {
      return;
    }
    initialize_buffer(A);
  }

  /**
   * Move constructor.
   * @param A matrix_cl to move
   */
  matrix_cl(matrix_cl<T>&& A)
      : buffer_cl_(std::move(A.buffer_cl_)),
        rows_(A.rows_),
        cols_(A.cols_),
        view_(A.view_),
        write_events_(std::move(A.write_events_)),
        read_events_(std::move(A.read_events_)) {}

  /**
   * Constructor for the matrix_cl that creates a copy of a std::vector of Eigen
   * matrices on the OpenCL device. Each matrix is flattened into one column
   * of the resulting matrix_cl. If a lvalue is passed to this constructor the
   * caller must make sure that the vector does not go out of scope before
   * copying is complete.
   *
   * That means `.wait()` must be called on the event associated on copying or
   * any other event that requires completion of this event. This can be done by
   * calling `.wait_for_write_events()` or `.wait_for_read_write_events()` on
   * this matrix or any matrix that is calculated from this one.
   *
   * @param A the vector of Eigen matrices
   *
   * @throw <code>std::invalid_argument</code> if the
   * matrices do not have matching dimensions
   * @throw <code>std::system_error</code> if the memory on the device could not
   * be allocated
   */
  template <typename Vec, require_std_vector_vt<is_eigen, Vec>...,
            require_st_same<Vec, T>...>
  explicit matrix_cl(Vec&& A) try : rows_(A.empty() ? 0 : A[0].size()),
                                    cols_(A.size()) {
    if (this->size() == 0) {
      return;
    }
    cl::Context& ctx = opencl_context.context();
    cl::CommandQueue& queue = opencl_context.queue();
    buffer_cl_ = cl::Buffer(ctx, CL_MEM_READ_WRITE, sizeof(T) * size());
    for (int i = 0, offset_size = 0; i < cols_; i++, offset_size += rows_) {
      check_size_match("matrix constructor", "input rows", A[i].size(),
                       "matrix_cl rows", rows_);
      cl::Event write_event;
      queue.enqueueWriteBuffer(
          buffer_cl_,
          opencl_context.in_order() || std::is_rvalue_reference<Vec&&>::value,
          sizeof(T) * offset_size, sizeof(T) * rows_, A[i].data(), nullptr,
          &write_event);
      this->add_write_event(write_event);
    }
  } catch (const cl::Error& e) {
    check_opencl_error("matrix constructor", e);
  }

  /**
   * Constructor for the matrix_cl that
   * only allocates the buffer on the OpenCL device.
   * Regardless of `partial_view`, whole matrix is stored.
   *
   * @param rows number of matrix rows, must be greater or equal to 0
   * @param cols number of matrix columns, must be greater or equal to 0
   * @param partial_view which part of the matrix is used
   *
   * @throw <code>std::system_error</code> if the memory on the device could not
   * be allocated
   *
   */
  matrix_cl(const int rows, const int cols,
            matrix_cl_view partial_view = matrix_cl_view::Entire)
      : rows_(rows), cols_(cols), view_(partial_view) {
    if (size() == 0) {
      return;
    }
    cl::Context& ctx = opencl_context.context();
    try {
      buffer_cl_
          = cl::Buffer(ctx, CL_MEM_READ_WRITE, sizeof(T) * rows_ * cols_);
    } catch (const cl::Error& e) {
      check_opencl_error("matrix constructor", e);
    }
  }

  /**
   * Constructor for the matrix_cl that creates a copy of the Eigen matrix or
   * Eigen expression on the OpenCL device. Regardless of `partial_view`, whole
   * matrix is stored. If a lvalue matrix is passed to this constructor the
   * caller must make sure that the matrix does not go out of scope before
   * copying is complete.
   *
   * That means `.wait()` must be called on the event associated on copying or
   * any other event that requires completion of this event. This can be done by
   * calling `.wait_for_write_events()` or `.wait_for_read_write_events()` on
   * this matrix or any matrix that is calculated from this one.
   *
   * @tparam Mat type of \c Eigen \c Matrix or expression
   * @param A the \c Eigen \c Matrix or expression
   * @param partial_view which part of the matrix is used
   *
   * @throw <code>std::system_error</code> if the memory on the device could not
   * be allocated
   */
  template <typename Mat, require_eigen_t<Mat>..., require_vt_same<Mat, T>...>
  explicit matrix_cl(Mat&& A,
                     matrix_cl_view partial_view = matrix_cl_view::Entire)
      : rows_(A.rows()), cols_(A.cols()), view_(partial_view) {
    using Mat_type = std::decay_t<ref_type_for_opencl_t<Mat>>;
    if (size() == 0) {
      return;
    }
    if (std::is_same<std::decay_t<Mat>, Mat_type>::value
        && (std::is_lvalue_reference<Mat>::value
            || is_eigen_contiguous_map<Mat>::value)) {
      // .eval)= is here just in case other branch is selected and A does not
      // have data directley accessible
      initialize_buffer(A.eval().data());
    } else {
      auto* A_heap = new Mat_type(std::move(A));
      try {
        cl::Event e = initialize_buffer(A_heap->data());
        // We set a callback that will delete the memory once copying is
        // complete. This event object does not hold the information about
        // callback. OpenCL implementation does. So nothing is lost as event
        // goes out of scope.
        e.setCallback(CL_COMPLETE, &delete_it<Mat_type>, A_heap);
      } catch (...) {
        delete A_heap;
        throw;
      }
    }
  }

  /**
   * Constructor for the matrix_cl that creates a copy of a scalar on the OpenCL
   * device. Regardless of `partial_view`, whole matrix is stored. If a lvalue
   * is passed to this constructor the caller must make sure that it does not go
   * out of scope before copying is complete.
   *
   * That means `.wait()` must be called on the event associated on copying or
   * any other event that requires completion of this event. This can be done by
   * calling `.wait_for_write_events()` or `.wait_for_read_write_events()` on
   * this matrix or any matrix that is calculated from this one.
   *
   * @param A the scalar
   * @param partial_view which part of the matrix is used
   *
   * @throw <code>std::system_error</code> if the memory on the device could not
   * be allocated
   */
  template <typename Scal,
            typename = require_same_t<T, std::remove_reference_t<Scal>>>
  explicit matrix_cl(Scal&& A,
                     matrix_cl_view partial_view = matrix_cl_view::Diagonal)
      : rows_(1), cols_(1), view_(partial_view) {
    initialize_buffer<std::is_rvalue_reference<Scal&&>::value>(&A);
  }

  /**
   * Construct a matrix_cl of size Nx1 from \c std::vector. If a lvalue is
   * passed to this constructor the caller must make sure that it does not go
   * out of scope before copying is complete.
   *
   * That means `.wait()` must be called on the event associated on copying or
   * any other event that requires completion of this event. This can be done by
   * calling `.wait_for_write_events()` or `.wait_for_read_write_events()` on
   * this matrix or any matrix that is calculated from this one.
   *
   * @param A Standard vector
   * @param partial_view which part of the matrix is used
   *
   * @throw <code>std::system_error</code> if the memory on the device could not
   * be allocated
   */
  template <typename Vec, require_std_vector_t<Vec>...,
            require_vt_same<Vec, T>...>
  explicit matrix_cl(Vec&& A,
                     matrix_cl_view partial_view = matrix_cl_view::Entire)
      : matrix_cl(std::forward<Vec>(A), A.size(), 1) {}

  /**
   * Construct from \c std::vector with given rows and columns. If a lvalue
   * is passed to this constructor the caller must make sure that it does not
   * go out of scope before copying is complete.
   *
   * That means `.wait()` must be called on the event associated on copying or
   * any other event that requires completion of this event. This can be done by
   * calling `.wait_for_write_events()` or `.wait_for_read_write_events()` on
   * this matrix or any matrix that is calculated from this one.
   *
   * @param A Standard vector
   * @param R Number of rows the matrix should have.
   * @param C Number of columns the matrix should have.
   * @param partial_view which part of the matrix is used
   *
   * @throw <code>std::system_error</code> if the memory on the device could not
   * be allocated
   */
  template <typename Vec, require_std_vector_t<Vec>...,
            require_vt_same<Vec, T>...>
  explicit matrix_cl(Vec&& A, const int& R, const int& C,
                     matrix_cl_view partial_view = matrix_cl_view::Entire)
      : rows_(R), cols_(C), view_(partial_view) {
    initialize_buffer_optionally_from_heap(std::forward<Vec>(A));
  }

  /**
   * Construct from \c array with given rows and columns. The caller
   * must make sure that data is not deleted before copying is complete.
   *
   * That means `.wait()` must be called on the event associated on copying or
   * any other event that requires completion of this event. This can be done by
   * calling `.wait_for_write_events()` or `.wait_for_read_write_events()` on
   * this matrix or any matrix that is calculated from this one.
   *
   * @param A array of doubles
   * @param R Number of rows the matrix should have.
   * @param C Number of columns the matrix should have.
   * @param partial_view which part of the matrix is used
   *
   * @throw <code>std::system_error</code> if the memory on the device could not
   * be allocated
   */
  template <typename U, require_same_t<T, U>...>
  explicit matrix_cl(const U* A, const int& R, const int& C,
                     matrix_cl_view partial_view = matrix_cl_view::Entire)
      : rows_(R), cols_(C), view_(partial_view) {
    initialize_buffer(A);
  }

  /**
   * Construct from a kernel generator expression. It evaluates the expression
   * into \c this.
   * @tparam Expr type of the expression
   * @param expression expression
   */
  template <typename Expr,
            require_all_kernel_expressions_and_none_scalar_t<Expr>* = nullptr>
  matrix_cl(const Expr& expression);  // NOLINT(runtime/explicit)

  /**
   * Move assignment operator.
   */
  matrix_cl<T>& operator=(matrix_cl<T>&& a) {
    view_ = a.view();
    rows_ = a.rows();
    cols_ = a.cols();
    this->wait_for_read_write_events();
    buffer_cl_ = std::move(a.buffer_cl_);
    write_events_ = std::move(a.write_events_);
    read_events_ = std::move(a.read_events_);
    return *this;
  }

  /**
   * Copy assignment operator.
   */
  matrix_cl<T>& operator=(const matrix_cl<T>& a) {
    this->view_ = a.view();
    this->rows_ = a.rows();
    this->cols_ = a.cols();
    if (a.size() == 0) {
      return *this;
    }
    this->wait_for_read_write_events();
    initialize_buffer(a);
    return *this;
  }

  /**
   * Assignment of a kernel generator expression evaluates the expression into
   * \c this.
   * @tparam Expr type of the expression
   * @param expression expression
   */
  template <typename Expr,
            require_all_kernel_expressions_and_none_scalar_t<Expr>* = nullptr>
  matrix_cl<T>& operator=(const Expr& expression);

 private:
  /**
   * Initializes the OpenCL buffer of this matrix by copying the data from given
   * buffer. Assumes that size of \c this is already set and matches the
   * buffer size. If \c in_order is false the caller must make sure that data
   * is not deleted before copying is complete.
   *
   * That means `.wait()` must be called on the event associated on copying or
   * any other event that requires completion of this event. This can be done by
   * calling `.wait_for_write_events()` or `.wait_for_read_write_events()` on
   * this matrix or any matrix that is calculated from this one.
   *
   * @tparam in_order whether copying must be done in order
   * @param A pointer to buffer
   * @return event for the copy
   */
  template <bool in_order = false>
  cl::Event initialize_buffer(const T* A) {
    cl::Event transfer_event;
    if (size() == 0) {
      return transfer_event;
    }
    cl::Context& ctx = opencl_context.context();
    cl::CommandQueue& queue = opencl_context.queue();
    try {
      buffer_cl_ = cl::Buffer(ctx, CL_MEM_READ_WRITE, sizeof(T) * size());
      queue.enqueueWriteBuffer(buffer_cl_,
                               opencl_context.in_order() || in_order, 0,
                               sizeof(T) * size(), A, nullptr, &transfer_event);
      this->add_write_event(transfer_event);
    } catch (const cl::Error& e) {
      check_opencl_error("initialize_buffer", e);
    }
    return transfer_event;
  }

  /**
   * Initializes the OpenCL buffer of this matrix by copying the data from given
   * object. Assumes that size of \c this is already set and matches the
   * buffer size. If the object is rvalue (temporary) it is first moved to heap
   * and callback is set to delete it after copying to OpenCL device is
   * complete. If a lvalue is passed to this function the caller must make
   * sure that input object does not go out of scope before copying is complete.
   *
   * That means `.wait()` must be called on the event associated on copying or
   * any other event that requires completion of this event. This can be done by
   * calling `.wait_for_write_events()` or `.wait_for_read_write_events()` on
   * this matrix or any matrix that is calculated from this one.
   *
   * @tparam U type of object
   * @param obj object
   * @return event for the copy
   */
  template <typename U>
  void initialize_buffer_optionally_from_heap(U&& obj) {
    using U_val = std::decay_t<U>;
    if (size() == 0) {
      return;
    }
    if (std::is_lvalue_reference<U>::value) {
      initialize_buffer(obj.data());
    } else {
      auto* obj_heap = new U_val(std::move(obj));
      try {
        cl::Event e = initialize_buffer(obj_heap->data());
        e.setCallback(CL_COMPLETE, &delete_it<U_val>, obj_heap);
      } catch (...) {
        delete obj_heap;
        throw;
      }
    }
  }

  /**
   * Initializes the OpenCL buffer of this matrix by copying the data from given
   * matrix_cl. Assumes that size of \c this is already set and matches the
   * size of given matrix.
   * @param A matrix_cl
   */
  void initialize_buffer(const matrix_cl<T>& A) {
    cl::Context& ctx = opencl_context.context();
    cl::CommandQueue queue = opencl_context.queue();
    try {
      buffer_cl_ = cl::Buffer(ctx, CL_MEM_READ_WRITE, sizeof(T) * this->size());
      cl::Event cstr_event;
      queue.enqueueCopyBuffer(A.buffer(), this->buffer(), 0, 0,
                              A.size() * sizeof(T), &A.write_events(),
                              &cstr_event);
      this->add_write_event(cstr_event);
      A.add_read_event(cstr_event);
    } catch (const cl::Error& e) {
      check_opencl_error("copy (OpenCL)->(OpenCL)", e);
    }
  }

  /**
   * Deletes the container. Used as a callback for OpenCL event.
   * @tparam U type of container
   * @param e cl_event handle
   * @param status status of event
   * @param container container to delete
   */
  template <typename U>
  static void delete_it(cl_event e, cl_int status, void* container) {
    delete static_cast<U*>(container);
  }
};

template <typename T>
using matrix_cl_prim = matrix_cl<T, require_arithmetic_t<T>>;

template <typename T>
using matrix_cl_fp = matrix_cl<T, require_floating_point_t<T>>;

/** @}*/

}  // namespace math
}  // namespace stan

#endif
#endif
