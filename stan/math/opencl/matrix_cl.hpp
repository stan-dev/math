#ifndef STAN_MATH_OPENCL_MATRIX_CL_HPP
#define STAN_MATH_OPENCL_MATRIX_CL_HPP
#ifdef STAN_OPENCL
#include <stan/math/opencl/opencl_context.hpp>
#include <stan/math/opencl/constants.hpp>
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
  /**
   * cl::Buffer provides functionality for working with the OpenCL buffer.
   * An OpenCL buffer allocates the memory in the device that
   * is provided by the context.
   */
  cl::Buffer buffer_cl_;
  const int rows_;
  const int cols_;
  mutable std::vector<cl::Event> write_events_;  // Tracks write jobs
  mutable std::vector<cl::Event> read_events_;   // Tracks reads

 public:
  typedef T type;
  // Forward declare the methods that work in place on the matrix
  template <TriangularViewCL triangular_view = TriangularViewCL::Entire>
  void zeros();
  template <TriangularMapCL triangular_map = TriangularMapCL::LowerToUpper>
  void triangular_transpose();
  template <TriangularViewCL triangular_view = TriangularViewCL::Entire>
  void sub_block(const matrix_cl<T, enable_if_arithmetic<T>>& A, size_t A_i,
                 size_t A_j, size_t this_i, size_t this_j, size_t nrows,
                 size_t ncols);
  int rows() const { return rows_; }

  int cols() const { return cols_; }

  int size() const { return rows_ * cols_; }

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

  matrix_cl(const matrix_cl<T>& A) : rows_(A.rows()), cols_(A.cols()) {
    if (A.size() == 0)
      return;
    cl::Context& ctx = opencl_context.context();
    cl::CommandQueue queue = opencl_context.queue();
    try {
      buffer_cl_ = cl::Buffer(ctx, CL_MEM_READ_WRITE, sizeof(T) * size());
      cl::Event cstr_event;
      queue.enqueueCopyBuffer(A.buffer(), this->buffer(), 0, 0,
                              A.size() * sizeof(T), &A.write_events(),
                              &cstr_event);
      this->add_write_event(cstr_event);
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
  matrix_cl(const int& rows, const int& cols) : rows_(rows), cols_(cols) {
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
   * Construct from Eigen Matrix
   *
   * @tparam R Rows of the @c Eigen @c Matrix
   * @tparam C Columns of the @c Eigen @c Matrix
   *
   * @param A The @c Eigen @c Matrix to move to OpenCL device.
   *
   * @throw <code>std::system_error</code> if the
   * matrices do not have matching dimensions
   */
  template <int R, int C>
  explicit matrix_cl(const Eigen::Matrix<T, R, C>& A)
      : rows_(A.rows()), cols_(A.cols()) {
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
   * Construct from @c std::vector with given rows and columns
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
   * Assign a @c matrix_cl to another
   */
  matrix_cl<T>& operator=(const matrix_cl<T>& a) {
    check_size_match("assignment of (OpenCL) matrices", "source.rows()",
                     a.rows(), "destination.rows()", rows());
    check_size_match("assignment of (OpenCL) matrices", "source.cols()",
                     a.cols(), "destination.cols()", cols());
    // Need to wait for all of matrices events before destroying old buffer
    this->wait_for_read_write_events();
    buffer_cl_ = a.buffer();
    return *this;
  }

  /**
   * Assign a @c matrix_cl of one arithmetic type to another
   */
  template <typename U, typename = enable_if_arithmetic<T>>
  matrix_cl<T>& operator=(const matrix_cl<U>& a) {
    check_size_match("assignment of (OpenCL) matrices", "source.rows()",
                     a.rows(), "destination.rows()", rows());
    check_size_match("assignment of (OpenCL) matrices", "source.cols()",
                     a.cols(), "destination.cols()", cols());
    // Need to wait for all of matrices events before destroying old buffer
    this->wait_for_read_write_events();
    buffer_cl_ = a.buffer();
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
