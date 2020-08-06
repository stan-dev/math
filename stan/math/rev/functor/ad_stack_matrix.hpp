#ifndef STAN_MATH_REV_FUNCTOR_AD_STACK_MATRIX_HPP
#define STAN_MATH_REV_FUNCTOR_AD_STACK_MATRIX_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/core.hpp>

namespace stan {
namespace math {

/**
 * Equivalent to `Eigen::Matrix`, except that the data is stored on AD stack.
 * That makes these objects triviali destructible and usable in `vari`s.
 *
 * @tparam MatrixType Eigen matrix type this works as (`MatrixXd`, `VectorXd`
 * ...)
 */
template <typename MatrixType>
class AD_stack_matrix : public Eigen::Map<MatrixType> {
 public:
  using Scalar = value_type_t<MatrixType>;
  static constexpr int RowsAtCompileTime = MatrixType::RowsAtCompileTime;
  static constexpr int ColsAtCompileTime = MatrixType::ColsAtCompileTime;

  /**
   * Default constructor.
   */
  AD_stack_matrix()
      : Eigen::Map<MatrixType>::Map(
            nullptr,
            RowsAtCompileTime == Eigen::Dynamic ? 0 : RowsAtCompileTime,
            ColsAtCompileTime == Eigen::Dynamic ? 0 : ColsAtCompileTime) {}

  /**
   * Constructs `AD_stack_matrix` with given number of rows and columns.
   * @param rows number of rows
   * @param cols number of columns
   */
  AD_stack_matrix(Eigen::Index rows, Eigen::Index cols)
      : Eigen::Map<MatrixType>::Map(
            ChainableStack::instance_->memalloc_.alloc_array<Scalar>(rows
                                                                     * cols),
            rows, cols) {}

  /**
   * Constructs `AD_stack_matrix` with given size. This only works if
   * `MatrixType` is row or col vector.
   * @param size number of elements
   */
  explicit AD_stack_matrix(Eigen::Index size)
      : Eigen::Map<MatrixType>::Map(
            ChainableStack::instance_->memalloc_.alloc_array<Scalar>(size),
            size) {}

  /**
   * Constructs `AD_stack_matrix` from an expression.
   * @param other expression
   */
  template <typename T, require_eigen_t<T>* = nullptr>
  AD_stack_matrix(const T& other)  // NOLINT
      : Eigen::Map<MatrixType>::Map(
            ChainableStack::instance_->memalloc_.alloc_array<Scalar>(
                other.size()),
            other.rows(), other.cols()) {
    *this = other;
  }

  /**
   * Copy constructor.
   * @param other matrix to copy from
   */
  AD_stack_matrix(const AD_stack_matrix<MatrixType>& other)
      : Eigen::Map<MatrixType>::Map(const_cast<Scalar*>(other.data()),
                                    other.rows(), other.cols()) {}

  using Eigen::Map<MatrixType>::operator=;

  /**
   * Copy assignment operator.
   * @param other matrix to copy from
   * @return `*this`
   */
  AD_stack_matrix& operator=(const AD_stack_matrix<MatrixType>& other) {
    // placement new changes what data map points to - there is no allocation
    new (this) Eigen::Map<MatrixType>(const_cast<Scalar*>(other.data()),
                                      other.rows(), other.cols());
    return *this;
  }

  /**
   * Assignment operator for assigning an expression.
   * @param other expression to evaluate into this
   * @return `*this`
   */
  template <typename T>
  AD_stack_matrix& operator=(const T& a) {
    // placement new changes what data map points to - there is no allocation
    new (this) Eigen::Map<MatrixType>(
        ChainableStack::instance_->memalloc_.alloc_array<Scalar>(a.size()),
        a.rows(), a.cols());
    Eigen::Map<MatrixType>::operator=(a);
    return *this;
  }

  /**
   * Resizes this object to given number of rows and columns.
   * @param rows number of rows
   * @param cols number of columns
   */
  void resize(Eigen::Index rows, Eigen::Index cols) {
    // placement new changes what data map points to - this new is not an
    // allocation
    new (this) Eigen::Map<MatrixType>(
        ChainableStack::instance_->memalloc_.alloc_array<Scalar>(rows * cols),
        rows, cols);
  }

  /**
   * Resizes this object to given size. This only works if `MatrixType` is row
   * or col vector.
   * @param size number of elements
   */
  void resize(Eigen::Index size) {
    // placement new changes what data map points to - this new is not an
    // allocation
    new (this) Eigen::Map<MatrixType>(
        ChainableStack::instance_->memalloc_.alloc_array<Scalar>(size), size);
  }
};

}  // namespace math
}  // namespace stan

#endif
