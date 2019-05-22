#ifndef STAN_MATH_OPENCL_MATRIX_CL_HPP
#define STAN_MATH_OPENCL_MATRIX_CL_HPP
#ifdef STAN_OPENCL
#include <stan/math/opencl/opencl_context.hpp>
#include <stan/math/opencl/constants.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/opencl/cache_copy.hpp>
#include <stan/math/opencl/err/check_opencl.hpp>
#include <stan/math/prim/scal/err/check_size_match.hpp>
#include <stan/math/prim/scal/err/domain_error.hpp>
#include <stan/math/prim/arr/fun/vec_concat.hpp>
#include <stan/math/rev/scal/fun/value_of_rec.hpp>
#include <CL/cl.hpp>
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
  mutable std::vector<cl::Event> write_events_;  // Tracks write jobs
  mutable std::vector<cl::Event> read_events_;   // Tracks reads

 public:
  // Forward declare the methods that work in place on the matrix
  template <TriangularViewCL triangular_view = TriangularViewCL::Entire>
  void zeros();
  template <TriangularMapCL triangular_map = TriangularMapCL::LowerToUpper>
  void triangular_transpose();
  template <TriangularViewCL triangular_view = TriangularViewCL::Entire>
  void sub_block(const matrix_cl& A, size_t A_i, size_t A_j, size_t this_i,
                 size_t this_j, size_t nrows, size_t ncols);
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

  const cl::Buffer& buffer() const { return oclBuffer_; }
  matrix_cl() : rows_(0), cols_(0) {}

  matrix_cl(cl::Buffer& A, const int R, const int C) : rows_(R), cols_(C) {
    oclBuffer_ = A;
  }
  template <typename T>
  explicit matrix_cl(std::vector<T> A) : rows_(1), cols_(A.size()) {
    if (A.size() == 0)
      return;
    // the context is needed to create the buffer object
    cl::Context& ctx = opencl_context.context();
    cl::CommandQueue& queue = opencl_context.queue();
    try {
      oclBuffer_
          = cl::Buffer(ctx, CL_MEM_READ_WRITE, sizeof(double) * A.size());
      queue.enqueueWriteBuffer(oclBuffer_, CL_TRUE, 0,
                               sizeof(double) * A.size(),
                               value_of_rec(A).data());
    } catch (const cl::Error& e) {
      check_opencl_error("matrix_cl(std::vector<T>) constructor", e);
    }
  }

  explicit matrix_cl(std::vector<double> A) : rows_(1), cols_(A.size()) {
    if (A.size() == 0)
      return;
    // the context is needed to create the buffer object
    cl::Context& ctx = opencl_context.context();
    cl::CommandQueue& queue = opencl_context.queue();
    try {
      oclBuffer_
          = cl::Buffer(ctx, CL_MEM_READ_WRITE, sizeof(double) * A.size());
      queue.enqueueWriteBuffer(oclBuffer_, CL_TRUE, 0,
                               sizeof(double) * A.size(), A.data());
    } catch (const cl::Error& e) {
      check_opencl_error("matrix_cl(std::vector<T>) constructor", e);
    }
  }

  template <typename T>
  explicit matrix_cl(std::vector<T> A, const int R, const int C)
      : rows_(R), cols_(C) {
    if (A.size() == 0)
      return;
    // the context is needed to create the buffer object
    cl::Context& ctx = opencl_context.context();
    cl::CommandQueue& queue = opencl_context.queue();
    try {
      oclBuffer_
          = cl::Buffer(ctx, CL_MEM_READ_WRITE, sizeof(double) * A.size());
      queue.enqueueWriteBuffer(oclBuffer_, CL_TRUE, 0,
                               sizeof(double) * A.size(),
                               value_of_rec(A).data());
    } catch (const cl::Error& e) {
      check_opencl_error("matrix_cl(std::vector<T>, rows, cols) constructor",
                         e);
    }
  }

  explicit matrix_cl(std::vector<double> A, const int R, const int C)
      : rows_(R), cols_(C) {
    if (A.size() == 0)
      return;
    // the context is needed to create the buffer object
    cl::Context& ctx = opencl_context.context();
    cl::CommandQueue& queue = opencl_context.queue();
    try {
      oclBuffer_
          = cl::Buffer(ctx, CL_MEM_READ_WRITE, sizeof(double) * A.size());
      queue.enqueueWriteBuffer(oclBuffer_, CL_TRUE, 0,
                               sizeof(double) * A.size(), A.data());
    } catch (const cl::Error& e) {
      check_opencl_error("matrix_cl(std::vector<T>, rows, cols) constructor",
                         e);
    }
  }

  matrix_cl(const matrix_cl& A) : rows_(A.rows()), cols_(A.cols()) {
    if (A.size() == 0)
      return;
    // the context is needed to create the buffer object
    cl::Context& ctx = opencl_context.context();
    cl::CommandQueue queue = opencl_context.queue();
    try {
      // creates a read&write object for "size" double values
      // in the provided context
      oclBuffer_ = cl::Buffer(ctx, CL_MEM_READ_WRITE, sizeof(double) * size());
      cl::Event cstr_event;
      queue.enqueueCopyBuffer(A.buffer(), this->buffer(), 0, 0,
                              A.size() * sizeof(double), &A.write_events(),
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
      oclBuffer_
          = cl::Buffer(ctx, CL_MEM_READ_WRITE, sizeof(double) * rows_ * cols_);
    } catch (const cl::Error& e) {
      check_opencl_error("matrix constructor", e);
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

  /**
   * Constructor for the matrix_cl that
   * creates a copy of a var type Eigen matrix on the GPU.
   *
   * @tparam R rows of matrix
   * @tparam C cols of matrix
   * @param A the Eigen matrix
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

  matrix_cl& operator=(const matrix_cl& a) {
    check_size_match("assignment of (OpenCL) matrices", "source.rows()",
                     a.rows(), "destination.rows()", rows());
    check_size_match("assignment of (OpenCL) matrices", "source.cols()",
                     a.cols(), "destination.cols()", cols());
    // Need to wait for all of matrices events before destroying old buffer
    this->wait_for_read_write_events();
    oclBuffer_ = a.buffer();
    return *this;
  }
};
/**
 * Represents a matrix of varis on the OpenCL device.
 *
 * The matrix data is stored in two separate matrix_cl
 * members for values (val_) and adjoints (adj_).
 */
template <TriangularViewCL triangular_view = TriangularViewCL::Entire>
class matrix_v_cl {
 private:
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
  matrix_v_cl(vari**& A, int M) : rows_(M), cols_(M), val_(M, M), adj_(M, M) {
    if (size() == 0)
      return;
    cl::Context& ctx = opencl_context.context();
    cl::CommandQueue& queue = opencl_context.queue();
    try {
      const int vari_size = M * (M + 1) / 2;
      std::vector<double> val_cpy(vari_size);
      std::vector<double> adj_cpy(vari_size);
      for (size_t j = 0; j < vari_size; ++j) {
        val_cpy[j] = A[j]->val_;
        adj_cpy[j] = A[j]->adj_;
      }
      matrix_cl packed_val(val_cpy);
      matrix_cl packed_adj(adj_cpy);
      queue.enqueueWriteBuffer(packed_val.buffer(), CL_TRUE, 0,
                               sizeof(double) * vari_size, val_cpy.data());
      queue.enqueueWriteBuffer(packed_adj.buffer(), CL_TRUE, 0,
                               sizeof(double) * vari_size, adj_cpy.data());
      stan::math::opencl_kernels::unpack(cl::NDRange(M, M), val_.buffer(),
                                         packed_val.buffer(), M, M,
                                         triangular_view);
      stan::math::opencl_kernels::unpack(cl::NDRange(M, M), adj_.buffer(),
                                         packed_adj.buffer(), M, M,
                                         triangular_view);
    } catch (const cl::Error& e) {
      check_opencl_error("matrix constructor", e);
    }
  }
};

// if the triangular view is entire dont unpack, just copy
template <>
inline matrix_v_cl<TriangularViewCL::Entire>::matrix_v_cl(vari**& A, int M)
    : rows_(M), cols_(M), val_(M, M), adj_(M, M) {
  if (size() == 0)
    return;
  const int vari_size = M * M;
  std::vector<double> val_cpy(vari_size);
  std::vector<double> adj_cpy(vari_size);
  for (size_t j = 0; j < vari_size; ++j) {
    val_cpy[j] = A[j]->val_;
    adj_cpy[j] = A[j]->adj_;
  }
  val_ = matrix_cl(val_cpy, M, M);
  val_ = matrix_cl(adj_cpy, M, M);
}

}  // namespace math
}  // namespace stan

#endif
#endif
