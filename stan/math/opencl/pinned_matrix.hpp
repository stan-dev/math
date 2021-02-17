#ifndef STAN_MATH_OPENCL_PINNED_MATRIX_HPP
#define STAN_MATH_OPENCL_PINNED_MATRIX_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/opencl/opencl_context.hpp>
#include <CL/cl2.hpp>

namespace stan {
namespace math {

/**
 * Equivalent to `Eigen::Matrix`, except that the data is stored in (hopefully
 * pinned) memory, allocated by OpenCL driver.
 *
 * If the OpenCL device uses the same memory as host, this memory can be used by
 * device without making copies. Otherwise copying to and from device using this
 * memory is usually faster than using non-pinned memory.
 *
 * @tparam MatrixType Eigen matrix type this works as (`MatrixXd`, `VectorXd`
 * ...)
 */
template <typename MatrixType>
class pinned_matrix : public Eigen::Map<MatrixType> {
 public:
  using Scalar = value_type_t<MatrixType>;
  using Base = Eigen::Map<MatrixType>;
  using PlainObject = std::decay_t<MatrixType>;
  static constexpr int RowsAtCompileTime = MatrixType::RowsAtCompileTime;
  static constexpr int ColsAtCompileTime = MatrixType::ColsAtCompileTime;

 protected:
  cl::Buffer buffer_;

  pinned_matrix(cl::Buffer&& b, Eigen::Index rows, Eigen::Index cols)
      : Base::Map(static_cast<Scalar*>(opencl_context.queue().enqueueMapBuffer(
                      b, true, CL_MAP_WRITE_INVALIDATE_REGION, 0,
                      sizeof(Scalar) * rows * cols)),
                  rows, cols),
        buffer_(std::move(b)) {}

  pinned_matrix(cl::Buffer&& b, Eigen::Index size)
      : Base::Map(static_cast<Scalar*>(opencl_context.queue().enqueueMapBuffer(
                      b, true, CL_MAP_WRITE_INVALIDATE_REGION, 0,
                      sizeof(Scalar) * size)),
                  size),
        buffer_(std::move(b)) {}

 public:
  /**
   * Default constructor.
   */
  pinned_matrix()
      : Base::Map(nullptr,
                  RowsAtCompileTime == Eigen::Dynamic ? 0 : RowsAtCompileTime,
                  ColsAtCompileTime == Eigen::Dynamic ? 0 : ColsAtCompileTime) {
  }

  /**
   * Constructs `pinned_matrix` with given number of rows and columns.
   * @param rows number of rows
   * @param cols number of columns
   */
  pinned_matrix(Eigen::Index rows, Eigen::Index cols)
      : pinned_matrix(
            cl::Buffer(opencl_context.context(), CL_MEM_ALLOC_HOST_PTR,
                       sizeof(Scalar) * rows * cols),
            rows, cols) {}

  /**
   * Constructs `pinned_matrix` with given size. This only works if
   * `MatrixType` is row or col vector.
   * @param size number of elements
   */
  explicit pinned_matrix(Eigen::Index size)
      : pinned_matrix(cl::Buffer(opencl_context.context(),
                                 CL_MEM_ALLOC_HOST_PTR, sizeof(Scalar) * size),
                      size) {}

  /**
   * Constructs `pinned_matrix` from an expression.
   * @param other expression
   */
  template <typename T, require_eigen_t<T>* = nullptr>
  pinned_matrix(const T& other)  // NOLINT
      : pinned_matrix(other.rows(), other.cols()) {
    Base::operator=(other);
  }

  /**
   * Copy constructor.
   * @param other matrix to copy from
   */
  pinned_matrix(const pinned_matrix<MatrixType>& other)
      : pinned_matrix(other.rows(), other.cols()) {
    Base::operator=(other);
  }

  /**
   * Move constructor.
   * @param other matrix to move from
   */
  pinned_matrix(pinned_matrix<MatrixType>&& other)
      : Base::Map(other.data(), other.rows(), other.cols()),
        buffer_(std::move(other.buffer_)) {}

  // without this using, compiler prefers combination of implicit construction
  // and copy assignment to the inherited operator when assigned an expression
  using Base::operator=;

  /**
   * Copy assignment operator.
   * @param other matrix to copy from
   * @return `*this`
   */
  pinned_matrix& operator=(const pinned_matrix<MatrixType>& other) {
    if (this->rows() != other.rows() || this->cols() != other.cols()) {
      buffer_ = cl::Buffer(opencl_context.context(), CL_MEM_ALLOC_HOST_PTR,
                           sizeof(Scalar) * other.size());
      // placement new changes what data map points to - only allocation happens
      // when creating new Buffer
      new (this)
          Base(static_cast<Scalar*>(opencl_context.queue().enqueueMapBuffer(
                   buffer_, true, CL_MAP_WRITE_INVALIDATE_REGION, 0,
                   sizeof(Scalar) * other.size())),
               other.rows(), other.cols());
    }
    Base::operator=(other);
    return *this;
  }

  /**
   * Move assignment operator.
   * @param other matrix to move from
   * @return `*this`
   */
  pinned_matrix& operator=(pinned_matrix<MatrixType>&& other) {
    // placement new changes what data map points to - there is no allocation
    new (this)
        Base(const_cast<Scalar*>(other.data()), other.rows(), other.cols());
    buffer_ = std::move(other.buffer_);
    return *this;
  }

  /**
   * Assignment operator for assigning an expression.
   * @param other expression to evaluate into this
   * @return `*this`
   */
  template <typename T>
  pinned_matrix& operator=(const T& other) {
    if (this->rows() != other.rows() || this->cols() != other.cols()) {
      buffer_ = cl::Buffer(opencl_context.context(), CL_MEM_ALLOC_HOST_PTR,
                           sizeof(Scalar) * other.size());
      // placement new changes what data map points to - only allocation happens
      // when creating new Buffer
      new (this)
          Base(static_cast<Scalar*>(opencl_context.queue().enqueueMapBuffer(
                   buffer_, true, CL_MAP_WRITE_INVALIDATE_REGION, 0,
                   sizeof(Scalar) * other.size())),
               other.rows(), other.cols());
    }
    Base::operator=(other);
    return *this;
  }
};

}  // namespace math
}  // namespace stan

#endif
#endif
