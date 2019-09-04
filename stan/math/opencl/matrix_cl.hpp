#ifndef STAN_MATH_OPENCL_MATRIX_CL_HPP
#define STAN_MATH_OPENCL_MATRIX_CL_HPP
#ifdef STAN_OPENCL
#include <stan/math/opencl/opencl_context.hpp>
#include <stan/math/opencl/matrix_cl_view.hpp>
#include <stan/math/opencl/err/check_opencl.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/arr/fun/vec_concat.hpp>
#include <stan/math/prim/scal/err/check_size_match.hpp>
#include <stan/math/prim/scal/err/domain_error.hpp>
#include <CL/cl.hpp>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

/**
 *  @file stan/math/opencl/matrix_cl.hpp
 *  @brief The matrix_cl class - allocates memory space on the OpenCL device,
 *    functions for transfering matrices to and from OpenCL devices
 */
namespace stan {
namespace math {

// Dummy class to instantiate matrix_cl to enable for specific types.
template <typename T, typename = void>
class matrix_cl {};

/**
 * Represents a matrix on the OpenCL device.
 *
 * @tparam T an arithmetic type for the type stored in the OpenCL buffer.
 */
template <typename T>
class matrix_cl<T, enable_if_arithmetic<T>> {
 private:
  cl::Buffer buffer_cl_;  // Holds the allocated memory on the device
  const int rows_;
  const int cols_;
  matrix_cl_view view_;  // Holds info on if matrix is a special type
  mutable std::vector<cl::Event> write_events_;  // Tracks write jobs
  mutable std::vector<cl::Event> read_events_;   // Tracks reads

 public:
  typedef T type;
  // Forward declare the methods that work in place on the matrix
  template <matrix_cl_view matrix_view = matrix_cl_view::Entire>
  void zeros();
  template <TriangularMapCL triangular_map = TriangularMapCL::LowerToUpper>
  void triangular_transpose();

  void sub_block(const matrix_cl<T, enable_if_arithmetic<T>>& A, size_t A_i,
                 size_t A_j, size_t this_i, size_t this_j, size_t nrows,
                 size_t ncols);
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
    cl::CommandQueue queue = opencl_context.queue();
    cl::Event copy_event;
    queue.enqueueBarrierWithWaitList(&this->write_events(), &copy_event);
    copy_event.wait();
    write_events_.clear();
    return;
  }

  /**
   * Waits for the read events and clears the read event stack.
   */
  inline void wait_for_read_events() const {
    cl::CommandQueue queue = opencl_context.queue();
    cl::Event copy_event;
    queue.enqueueBarrierWithWaitList(&this->read_events(), &copy_event);
    copy_event.wait();
    read_events_.clear();
    return;
  }

  /**
   * Waits for read and write events to finish and clears the read, write, and
   * read/write event stacks.
   */
  inline void wait_for_read_write_events() const {
    cl::CommandQueue queue = opencl_context.queue();
    cl::Event copy_event;
    const std::vector<cl::Event> mat_events = this->read_write_events();
    queue.enqueueBarrierWithWaitList(&mat_events, &copy_event);
    copy_event.wait();
    read_events_.clear();
    write_events_.clear();
    return;
  }

  const cl::Buffer& buffer() const { return buffer_cl_; }
  cl::Buffer& buffer() { return buffer_cl_; }
  matrix_cl() : rows_(0), cols_(0) {}

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

  matrix_cl(const matrix_cl<T>& A)
      : rows_(A.rows()), cols_(A.cols()), view_(A.view()) {
    if (A.size() == 0)
      return;
    this->wait_for_read_write_events();
    cl::Context& ctx = opencl_context.context();
    cl::CommandQueue queue = opencl_context.queue();
    try {
      buffer_cl_ = cl::Buffer(ctx, CL_MEM_READ_WRITE, sizeof(T) * size());
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

  matrix_cl(matrix_cl<T>&& A)
      : buffer_cl_(A.buffer_cl_),
        rows_(A.rows_),
        cols_(A.cols_),
        view_(A.view_),
        write_events_(std::move(A.write_events_)),
        read_events_(std::move(A.read_events_)) {}

  /**
   * Constructor for the matrix_cl that
   * creates a copy of the Eigen matrix on the OpenCL device.
   *
   *
   * @tparam R row type
   * @tparam C column type
   * @param A the Eigen matrix
   *
   * @throw <code>std::invalid_argument</code> if the
   * matrices do not have matching dimensions
   */
  template <int R, int C>
  explicit matrix_cl(const std::vector<Eigen::Matrix<T, R, C>>& A) try
      : rows_(A.empty() ? 0 : A[0].size()),
        cols_(A.size()) {
    if (this->size() == 0)
      return;
    cl::Context& ctx = opencl_context.context();
    cl::CommandQueue& queue = opencl_context.queue();
    // creates the OpenCL buffer to copy the Eigen
    // matrix to the OpenCL device
    buffer_cl_ = cl::Buffer(ctx, CL_MEM_READ_WRITE, sizeof(T) * size());
    for (int i = 0, offset_size = 0; i < cols_; i++, offset_size += rows_) {
      check_size_match("matrix constructor", "input rows", A[i].size(),
                       "matrix_cl rows", rows_);
      /**
       * Writes the contents of A[i] to the OpenCL buffer
       * starting at the offset sizeof(double)*start.
       * CL_TRUE denotes that the call is blocking as
       * we do not want to execute any further kernels
       * on the device until we are sure that the data
       * is finished transfering
       */
      cl::Event write_event;
      queue.enqueueWriteBuffer(
          buffer_cl_, CL_FALSE, sizeof(double) * offset_size,
          sizeof(double) * rows_, A[i].data(), NULL, &write_event);
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
   * @throw <code>std::system_error</code> if the
   * matrices do not have matching dimensions
   *
   */
  matrix_cl(const int& rows, const int& cols,
            matrix_cl_view partial_view = matrix_cl_view::Entire)
      : rows_(rows), cols_(cols), view_(partial_view) {
    if (size() == 0) {
      return;
    }
    cl::Context& ctx = opencl_context.context();
    try {
      // creates the OpenCL buffer of the provided size
      buffer_cl_
          = cl::Buffer(ctx, CL_MEM_READ_WRITE, sizeof(T) * rows_ * cols_);
    } catch (const cl::Error& e) {
      check_opencl_error("matrix constructor", e);
    }
  }

  /**
   * Constructor for the matrix_cl that
   * creates a copy of the Eigen matrix on the OpenCL device.
   * Regardless of `partial_view`, whole matrix is stored.
   *
   * @tparam T type of data in the \c Eigen \c Matrix
   * @param A the \c Eigen \c Matrix
   * @param partial_view which part of the matrix is used
   *
   * @throw <code>std::system_error</code> if the
   * matrices do not have matching dimensions
   */
  explicit matrix_cl(
      const Eigen::Ref<const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>&
          A,
      matrix_cl_view partial_view = matrix_cl_view::Entire)
      : rows_(A.rows()), cols_(A.cols()), view_(partial_view) {
    if (size() == 0) {
      return;
    }
    cl::Context& ctx = opencl_context.context();
    cl::CommandQueue& queue = opencl_context.queue();
    try {
      buffer_cl_ = cl::Buffer(ctx, CL_MEM_READ_WRITE, sizeof(T) * A.size());
      cl::Event transfer_event;
      queue.enqueueWriteBuffer(buffer_cl_, CL_FALSE, 0, sizeof(T) * A.size(),
                               A.data(), NULL, &transfer_event);
      this->add_write_event(transfer_event);
    } catch (const cl::Error& e) {
      check_opencl_error("matrix constructor", e);
    }
  }

  /**
   * Construct from \c std::vector with given rows and columns
   *
   * @param A Standard vector
   * @param R Number of rows the matrix should have.
   * @param C Number of columns the matrix should have.
   * @throw <code>std::system_error</code> if the
   * matrices do not have matching dimensions
   */
  explicit matrix_cl(const std::vector<T>& A, const int& R, const int& C)
      : rows_(R), cols_(C) {
    if (size() == 0) {
      return;
    }
    cl::Context& ctx = opencl_context.context();
    cl::CommandQueue& queue = opencl_context.queue();
    try {
      buffer_cl_ = cl::Buffer(ctx, CL_MEM_READ_WRITE, sizeof(T) * A.size());
      cl::Event transfer_event;
      queue.enqueueWriteBuffer(buffer_cl_, CL_FALSE, 0, sizeof(T) * A.size(),
                               A.data(), NULL, &transfer_event);
      this->add_write_event(transfer_event);
    } catch (const cl::Error& e) {
      check_opencl_error("matrix constructor", e);
    }
  }

  /**
   * Constructs a const matrix_cl that contains a copy of the Eigen matrix on
   * the OpenCL device. If the matrix already has a cached copy on the device,
   * the cache is used and no copying is done. Changing the resulting matrix_cl
   * would change cache, so do not do it! If changes are needed a copy must be
   * made.
   *
   * @tparam R row type of input matrix
   * @tparam C column type of input matrix
   * @param A the Eigen matrix
   * @param partial_view which part of the matrix is used
   */
  template <int R, int C>
  static matrix_cl<T> constant(const Eigen::Matrix<T, R, C>& A,
                               matrix_cl_view partial_view
                               = matrix_cl_view::Entire) {
#ifndef STAN_OPENCL_NOCACHE
    if (A.opencl_buffer_() != NULL) {
      return matrix_cl<T>(A.opencl_buffer_, A.rows(), A.cols(), partial_view);
    } else {
      matrix_cl<T> res(A, partial_view);
      A.opencl_buffer_ = res.buffer();
      return res;
    }
#else
    return matrix_cl<T>(A, partial_view);
#endif
  }

  /**
   * Constructs a const matrix_cl that contains a single value on the OpenCL
   * device.
   *
   * @param A the value
   * @param partial_view which part of the matrix is used
   */
  static matrix_cl<T> constant(T A, matrix_cl_view partial_view
                                    = matrix_cl_view::Entire) {
    return matrix_cl<T>(A);
  }

  /**
   * Construct from \c array of doubles with given rows and columns
   *
   * @param A array of doubles
   * @param R Number of rows the matrix should have.
   * @param C Number of columns the matrix should have.
   * @param partial_view which part of the matrix is used
   * @throw <code>std::system_error</code> if the
   * matrices do not have matching dimensions
   */
  explicit matrix_cl(const double* A, const int& R, const int& C,
                     matrix_cl_view partial_view = matrix_cl_view::Entire)
      : rows_(R), cols_(C), view_(partial_view) {
    if (size() == 0) {
      return;
    }
    cl::Context& ctx = opencl_context.context();
    cl::CommandQueue& queue = opencl_context.queue();
    try {
      buffer_cl_ = cl::Buffer(ctx, CL_MEM_READ_WRITE, sizeof(T) * size());
      cl::Event transfer_event;
      queue.enqueueWriteBuffer(buffer_cl_, CL_FALSE, 0, sizeof(T) * size(), A,
                               NULL, &transfer_event);
      this->add_write_event(transfer_event);
    } catch (const cl::Error& e) {
      check_opencl_error("matrix constructor", e);
    }
  }

  /**
   * Assign a \c matrix_cl to another
   */
  matrix_cl<T>& operator=(const matrix_cl<T>& a) {
    check_size_match("assignment of (OpenCL) matrices", "source.rows()",
                     a.rows(), "destination.rows()", rows());
    check_size_match("assignment of (OpenCL) matrices", "source.cols()",
                     a.cols(), "destination.cols()", cols());
    if (a.size() == 0)
      return *this;
    view_ = a.view();
    this->wait_for_read_write_events();
    cl::CommandQueue queue = opencl_context.queue();
    try {
      cl::Event copy_event;
      queue.enqueueCopyBuffer(a.buffer(), this->buffer(), 0, 0,
                              a.size() * sizeof(T), &a.write_events(),
                              &copy_event);
      this->add_write_event(copy_event);
      a.add_read_event(copy_event);
    } catch (const cl::Error& e) {
      check_opencl_error("copy (OpenCL)->(OpenCL)", e);
    }
    return *this;
  }

  /**
   * Move a \c matrix_cl to another
   */
  matrix_cl<T>& operator=(matrix_cl<T>&& a) {
    check_size_match("move of (OpenCL) matrix", "source.rows()", a.rows(),
                     "destination.rows()", rows());
    check_size_match("move of (OpenCL) matrix", "source.cols()", a.cols(),
                     "destination.cols()", cols());
    // Need to wait for all of matrices events before destroying old buffer
    this->wait_for_read_write_events();
    buffer_cl_ = a.buffer();
    view_ = a.view();
    write_events_ = std::move(a.write_events_);
    read_events_ = std::move(a.read_events_);
    return *this;
  }

  /**
   * Assign a \c matrix_cl of one arithmetic type to another
   */
  template <typename U, typename = enable_if_arithmetic<U>>
  matrix_cl<T>& operator=(const matrix_cl<U>& a) {
    check_size_match("assignment of (OpenCL) matrices", "source.rows()",
                     a.rows(), "destination.rows()", rows());
    check_size_match("assignment of (OpenCL) matrices", "source.cols()",
                     a.cols(), "destination.cols()", cols());
    // Need to wait for all of matrices events before destroying old buffer
    this->wait_for_read_write_events();
    buffer_cl_ = a.buffer();
    view_ = a.view();
    return *this;
  }
};

template <typename T>
using matrix_cl_prim = matrix_cl<T, enable_if_arithmetic<T>>;

template <typename T>
using matrix_cl_fp = matrix_cl<T, enable_if_floating_point<T>>;

}  // namespace math
}  // namespace stan

#endif
#endif
