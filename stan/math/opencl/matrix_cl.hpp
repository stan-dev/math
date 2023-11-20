#ifndef STAN_MATH_OPENCL_MATRIX_CL_HPP
#define STAN_MATH_OPENCL_MATRIX_CL_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/prim/size.hpp>
#include <stan/math/opencl/err/check_opencl.hpp>
#include <stan/math/prim/err/check_size_match.hpp>
#include <stan/math/opencl/opencl_context.hpp>
#include <stan/math/opencl/ref_type_for_opencl.hpp>
#include <stan/math/opencl/matrix_cl_view.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/vec_concat.hpp>
#include <CL/opencl.hpp>
#include <tbb/concurrent_vector.h>
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

// forward declare
template <typename T>
class arena_matrix_cl;

template <typename>
class matrix_cl;

/**
 * Represents an arithmetic matrix on the OpenCL device.
 * @tparam T an arithmetic type for the type stored in the OpenCL buffer.
 */
template <typename T>
class matrix_cl : public matrix_cl_base {
 private:
  cl::Buffer buffer_cl_;  // Holds the allocated memory on the device
  int rows_{0};           // Number of rows.
  int cols_{0};           // Number of columns.
  // Holds info on if matrix is a special type
  matrix_cl_view view_{matrix_cl_view::Entire};
  mutable tbb::concurrent_vector<cl::Event> write_events_;  // Tracks write jobs
  mutable tbb::concurrent_vector<cl::Event> read_events_;   // Tracks reads

 public:
  using Scalar = T;  // Underlying type of the matrix
  using type = T;    // Underlying type of the matrix
  // Forward declare the methods that work in place on the matrix
  template <matrix_cl_view matrix_view = matrix_cl_view::Entire>
  inline void zeros_strict_tri();

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
  inline const tbb::concurrent_vector<cl::Event>& write_events() const {
    return write_events_;
  }

  /**
   * Get the events from the event stacks.
   * @return The read/write event stack.
   */
  inline const tbb::concurrent_vector<cl::Event>& read_events() const {
    return read_events_;
  }

  /**
   * Get the events from the event stacks.
   * @return The read/write event stack.
   */
  inline const tbb::concurrent_vector<cl::Event> read_write_events() const {
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
  matrix_cl(const cl::Buffer& A, const int R, const int C,
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
    buffer_cl_ = cl::Buffer(opencl_context.context(), CL_MEM_READ_WRITE,
                            sizeof(T) * this->size());
    initialize_buffer_cl(A);
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
   * Constructor from `arena_matrix_cl`.
   * @param A matrix_cl to move
   */
  // defined in rev/arena_matrix_cl.hpp
  matrix_cl(const arena_matrix_cl<T>& A);  // NOLINT(runtime/explicit)

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
  template <typename Vec, require_std_vector_vt<is_eigen, Vec>* = nullptr,
            require_st_same<Vec, T>* = nullptr>
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
    if (this->size() == 0) {
      return;
    }
    cl::Context& ctx = opencl_context.context();
    try {
      int flags = CL_MEM_READ_WRITE;
      if (opencl_context.device()[0].getInfo<CL_DEVICE_HOST_UNIFIED_MEMORY>()) {
        flags |= CL_MEM_ALLOC_HOST_PTR;
      }
      buffer_cl_ = cl::Buffer(ctx, flags, sizeof(T) * rows_ * cols_);
    } catch (const cl::Error& e) {
      check_opencl_error("matrix constructor", e);
    }
  }

  /**
   * Constructor for the matrix_cl that creates a copy of the Eigen matrix or
   * Eigen expression on the OpenCL device. Regardless of `partial_view`, whole
   * matrix is stored.
   *
   * If a lvalue matrix or a map is passed to this constructor, it might be
   * directly used by the device. The caller must make sure that the matrix (map
   * data) does not go out of scope as long as this `matrix_cl` is in use
   * (`std::move`-ing it or using raw `buffer()` also counts as in use).
   *
   * @tparam Mat type of \c Eigen \c Matrix or expression
   * @param A the \c Eigen \c Matrix or expression
   * @param partial_view which part of the matrix is used
   *
   * @throw <code>std::system_error</code> if the memory on the device could not
   * be allocated
   */
  template <typename Mat, require_eigen_t<Mat>* = nullptr,
            require_vt_same<Mat, T>* = nullptr>
  explicit matrix_cl(Mat&& A,
                     matrix_cl_view partial_view = matrix_cl_view::Entire)
      : rows_(A.rows()), cols_(A.cols()), view_(partial_view) {
    using Mat_type = std::decay_t<ref_type_for_opencl_t<Mat>>;
    if (this->size() == 0) {
      return;
    }
    initialize_buffer_no_heap_if<
        std::is_same<std::decay_t<Mat>, Mat_type>::value
        && (std::is_lvalue_reference<Mat>::value
            || is_eigen_contiguous_map<Mat>::value)>(A);
  }

  /**
   * Constructor for the matrix_cl that creates a copy of a scalar on the OpenCL
   * device. Regardless of `partial_view`, whole matrix is stored.
   *
   * If a lvalue is passed to this constructor, it might be directly used by the
   * device. The caller must make sure that it does not go out of scope as long
   * as this `matrix_cl` is in use
   * (`std::move`-ing it or using raw `buffer()` also counts as in use).
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
    initialize_buffer<std::is_rvalue_reference<Scal&&>::value>(
        const_cast<const std::decay_t<Scal>*>(&A));
  }

  /**
   * Construct a matrix_cl of size Nx1 from \c std::vector.
   *
   * If a lvalue is passed to this constructor, it might be directly used by the
   * device. The caller must make sure that it does not go out of scope as long
   * as this `matrix_cl` is in use
   * (`std::move`-ing it or using raw `buffer()` also counts as in use).
   *
   * @param A Standard vector
   * @param partial_view which part of the matrix is used
   *
   * @throw <code>std::system_error</code> if the memory on the device could not
   * be allocated
   */
  template <typename Vec, require_std_vector_t<Vec>* = nullptr,
            require_vt_same<Vec, T>* = nullptr>
  explicit matrix_cl(Vec&& A,
                     matrix_cl_view partial_view = matrix_cl_view::Entire)
      : matrix_cl(std::forward<Vec>(A), A.size(), 1) {}

  /**
   * Construct from \c std::vector with given rows and columns.
   *
   * If a lvalue is passed to this constructor, it might be directly used by the
   * device. The caller must make sure that it does not go out of scope as long
   * as this `matrix_cl` is in use `std::move`-ing it or using raw `buffer()`
   * also counts as in use).
   *
   * @param A Standard vector
   * @param R Number of rows the matrix should have.
   * @param C Number of columns the matrix should have.
   * @param partial_view which part of the matrix is used
   *
   * @throw <code>std::system_error</code> if the memory on the device could not
   * be allocated
   */
  template <typename Vec, require_std_vector_t<Vec>* = nullptr,
            require_vt_same<Vec, T>* = nullptr>
  explicit matrix_cl(Vec&& A, const int& R, const int& C,
                     matrix_cl_view partial_view = matrix_cl_view::Entire)
      : rows_(R), cols_(C), view_(partial_view) {
    initialize_buffer_no_heap_if<std::is_lvalue_reference<Vec>::value>(A);
  }

  /**
   * Construct from \c array with given rows and columns.
   *
   * The memory might be directly used by the device. The caller must make sure
   * that it does not go out of scope as long as this `matrix_cl` is in use
   * (`std::move`-ing it or using raw `buffer()` also counts as in use).
   *
   * @param A array of doubles
   * @param R Number of rows the matrix should have.
   * @param C Number of columns the matrix should have.
   * @param partial_view which part of the matrix is used
   *
   * @throw <code>std::system_error</code> if the memory on the device could not
   * be allocated
   */
  template <typename U, require_same_t<T, U>* = nullptr>
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
  // defined in kernel_generator/matrix_cl_conversion.hpp
  template <typename Expr,
            require_all_kernel_expressions_and_none_scalar_t<Expr>* = nullptr,
            require_not_matrix_cl_t<Expr>* = nullptr>
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
    if (a.size() == 0) {
      this->rows_ = a.rows();
      this->cols_ = a.cols();
      return *this;
    }
    this->wait_for_read_write_events();
    if (this->size() != a.size()) {
      buffer_cl_ = cl::Buffer(opencl_context.context(), CL_MEM_READ_WRITE,
                              sizeof(T) * a.size());
    }
    this->rows_ = a.rows();
    this->cols_ = a.cols();
    initialize_buffer_cl(a);
    return *this;
  }

  /**
   * Assignment of a kernel generator expression evaluates the expression into
   * \c this.
   * @tparam Expr type of the expression
   * @param expression expression
   */
  // defined in kernel_generator/matrix_cl_conversion.hpp
  template <typename Expr,
            require_all_kernel_expressions_and_none_scalar_t<Expr>* = nullptr,
            require_not_matrix_cl_t<Expr>* = nullptr>
  matrix_cl<T>& operator=(const Expr& expression);

  /**
   * Assignment of `arena_matrix_cl<T>`.
   * @tparam T type of matrix
   * @param other matrix
   */
  // defined in rev/arena_matrix_cl.hpp
  matrix_cl<T>& operator=(const arena_matrix_cl<T>& other);

  /**
   * Evaluates `this`. This is a no-op.
   * @return `*this`
   */
  const matrix_cl<T>& eval() const& { return *this; }
  matrix_cl<T> eval() && { return std::move(*this); }

  /**
   * Destructor waits for write events to prevent any kernels from writing
   * memory that has already been reused.
   */
  ~matrix_cl() { wait_for_read_write_events(); }

 private:
  /**
   * Initializes the OpenCL buffer of this matrix by copying the data from given
   * buffer. Assumes that size of \c this is already set and matches the
   * buffer size.
   *
   * The caller must make sure that data is not deleted as long as
   * this `matrix_cl` is in use (`std::move`-ing it or using raw `buffer()` also
   * counts as in use).
   *
   * @tparam in_order whether copying must be done in order
   * efficiently use it directly
   * @param A pointer to buffer
   * @return event for the copy
   */
  template <bool in_order = false>
  cl::Event initialize_buffer(const T* A) {
    cl::Event transfer_event;
    if (this->size() == 0) {
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

  template <bool in_order = false>
  cl::Event initialize_buffer(T* A) {
    cl::Event transfer_event;
    if (this->size() == 0) {
      return transfer_event;
    }
    cl::Context& ctx = opencl_context.context();
    cl::CommandQueue& queue = opencl_context.queue();
    try {
      if (opencl_context.device()[0].getInfo<CL_DEVICE_HOST_UNIFIED_MEMORY>()) {
        constexpr auto copy_or_share
            = CL_MEM_COPY_HOST_PTR * INTEGRATED_OPENCL
              | (CL_MEM_USE_HOST_PTR * !INTEGRATED_OPENCL);
        buffer_cl_
            = cl::Buffer(ctx, CL_MEM_READ_WRITE | copy_or_share,
                         sizeof(T) * size(), A);  // this is always synchronous
      } else {
        buffer_cl_ = cl::Buffer(ctx, CL_MEM_READ_WRITE, sizeof(T) * size());
        queue.enqueueWriteBuffer(
            buffer_cl_, opencl_context.in_order() || in_order, 0,
            sizeof(T) * size(), A, nullptr, &transfer_event);
        this->add_write_event(transfer_event);
      }
    } catch (const cl::Error& e) {
      check_opencl_error("initialize_buffer", e);
    }
    return transfer_event;
  }

  /**
   * Initializes the OpenCL buffer of this matrix by copying the data from given
   * object. Assumes that size of \c this is already set and matches the
   * buffer size. If No_heap is false the object is first moved to heap
   * and callback is set to delete it after copying to OpenCL device is
   * complete. Otherwise the caller must make sure that input object is not
   * deleted as long as this `matrix_cl` is in use (`std::move`-ing it or using
   * raw `buffer()` also counts).
   *
   * @tparam No_heap whether to move the object to heap first
   * @tparam U type of object
   * @param obj object
   * @return event for the copy
   */
  template <bool No_heap, typename U, std::enable_if_t<No_heap>* = nullptr>
  void initialize_buffer_no_heap_if(U&& obj) {
    if (this->size() == 0) {
      return;
    }
    initialize_buffer(obj.data());
  }
  // we need separate overloads as obj.data() might not be available when second
  // overload is called.
  template <bool No_heap, typename U, std::enable_if_t<!No_heap>* = nullptr>
  void initialize_buffer_no_heap_if(U&& obj) {
    using U_val = std::decay_t<ref_type_for_opencl_t<U>>;
    if (this->size() == 0) {
      return;
    }
    auto* obj_heap = new U_val(std::move(obj));
    try {
      cl::Event e = initialize_buffer(obj_heap->data());
      if (opencl_context.device()[0].getInfo<CL_DEVICE_HOST_UNIFIED_MEMORY>()) {
        buffer_cl_.setDestructorCallback(&delete_it_destructor<U_val>,
                                         obj_heap);
      } else {
        e.setCallback(CL_COMPLETE, &delete_it_event<U_val>, obj_heap);
      }
    } catch (...) {
      delete obj_heap;
      throw;
    }
  }

  /**
   * Initializes the OpenCL buffer of this matrix by copying the data from given
   * matrix_cl. Assumes that size of \c this is already set and matches the
   * size of given matrix.
   * @param A matrix_cl
   */
  void initialize_buffer_cl(const matrix_cl<T>& A) {
    cl::Event cstr_event;
    std::vector<cl::Event>* dep_events = new std::vector<cl::Event>(
        A.write_events().begin(), A.write_events().end());
    try {
      opencl_context.queue().enqueueCopyBuffer(A.buffer(), this->buffer(), 0, 0,
                                               A.size() * sizeof(T), dep_events,
                                               &cstr_event);
      if (opencl_context.device()[0].getInfo<CL_DEVICE_HOST_UNIFIED_MEMORY>()) {
        buffer_cl_.setDestructorCallback(
            &delete_it_destructor<std::vector<cl::Event>>, dep_events);
      } else {
        cstr_event.setCallback(
            CL_COMPLETE, &delete_it_event<std::vector<cl::Event>>, dep_events);
      }
      this->add_write_event(cstr_event);
      A.add_read_event(cstr_event);
    } catch (const cl::Error& e) {
      delete dep_events;
      check_opencl_error("copy (OpenCL)->(OpenCL)", e);
    } catch (...) {
      delete dep_events;
      throw;
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
  static void delete_it_event(cl_event e, cl_int status, void* container) {
    delete static_cast<U*>(container);
  }

  /**
   * Deletes the container. Used as a callback for destruction of `cl::Buffer`.
   * @tparam U type of container
   * @param buff buffer that is being destructed
   * @param container container to delete
   */
  template <typename U>
  static void delete_it_destructor(cl_mem buff, void* container) {
    delete static_cast<U*>(container);
  }
};
/** @}*/

}  // namespace math
}  // namespace stan

#endif
#endif
