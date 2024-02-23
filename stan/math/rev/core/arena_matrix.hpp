#ifndef STAN_MATH_REV_CORE_ARENA_MATRIX_HPP
#define STAN_MATH_REV_CORE_ARENA_MATRIX_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/core/chainable_alloc.hpp>
#include <stan/math/rev/core/chainablestack.hpp>

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
class arena_matrix : public Eigen::Map<MatrixType> {
 public:
  using Scalar = value_type_t<MatrixType>;
  using Base = Eigen::Map<MatrixType>;
  using PlainObject = std::decay_t<MatrixType>;
  static constexpr int RowsAtCompileTime = MatrixType::RowsAtCompileTime;
  static constexpr int ColsAtCompileTime = MatrixType::ColsAtCompileTime;

  /**
   * Default constructor.
   */
  arena_matrix()
      : Base::Map(nullptr,
                  RowsAtCompileTime == Eigen::Dynamic ? 0 : RowsAtCompileTime,
                  ColsAtCompileTime == Eigen::Dynamic ? 0 : ColsAtCompileTime) {
  }

  /**
   * Constructs `arena_matrix` with given number of rows and columns.
   * @param rows number of rows
   * @param cols number of columns
   */
  arena_matrix(Eigen::Index rows, Eigen::Index cols)
      : Base::Map(
          ChainableStack::instance_->memalloc_.alloc_array<Scalar>(rows * cols),
          rows, cols) {}

  /**
   * Constructs `arena_matrix` with given size. This only works if
   * `MatrixType` is row or col vector.
   * @param size number of elements
   */
  explicit arena_matrix(Eigen::Index size)
      : Base::Map(
          ChainableStack::instance_->memalloc_.alloc_array<Scalar>(size),
          size) {}

  /**
   * Constructs `arena_matrix` from an expression.
   * @param other expression
   */
  template <typename T, require_eigen_t<T>* = nullptr>
  arena_matrix(const T& other)  // NOLINT
      : Base::Map(
          ChainableStack::instance_->memalloc_.alloc_array<Scalar>(
              other.size()),
          (RowsAtCompileTime == 1 && T::ColsAtCompileTime == 1)
                  || (ColsAtCompileTime == 1 && T::RowsAtCompileTime == 1)
              ? other.cols()
              : other.rows(),
          (RowsAtCompileTime == 1 && T::ColsAtCompileTime == 1)
                  || (ColsAtCompileTime == 1 && T::RowsAtCompileTime == 1)
              ? other.rows()
              : other.cols()) {
    *this = other;
  }

  /**
   * Constructs `arena_matrix` from an expression. This makes an assumption that
   * any other `Eigen::Map` also contains memory allocated in the arena.
   * @param other expression
   */
  arena_matrix(const Base& other)  // NOLINT
      : Base::Map(other) {}

  /**
   * Copy constructor.
   * @param other matrix to copy from
   */
  arena_matrix(const arena_matrix<MatrixType>& other)
      : Base::Map(const_cast<Scalar*>(other.data()), other.rows(),
                  other.cols()) {}

  // without this using, compiler prefers combination of implicit construction
  // and copy assignment to the inherited operator when assigned an expression
  using Base::operator=;

  /**
   * Copy assignment operator.
   * @param other matrix to copy from
   * @return `*this`
   */
  arena_matrix& operator=(const arena_matrix<MatrixType>& other) {
    // placement new changes what data map points to - there is no allocation
    new (this)
        Base(const_cast<Scalar*>(other.data()), other.rows(), other.cols());
    return *this;
  }

  /**
   * Assignment operator for assigning an expression.
   * @param a expression to evaluate into this
   * @return `*this`
   */
  template <typename T>
  arena_matrix& operator=(const T& a) {
    // do we need to transpose?
    if ((RowsAtCompileTime == 1 && T::ColsAtCompileTime == 1)
        || (ColsAtCompileTime == 1 && T::RowsAtCompileTime == 1)) {
      // placement new changes what data map points to - there is no allocation
      new (this) Base(
          ChainableStack::instance_->memalloc_.alloc_array<Scalar>(a.size()),
          a.cols(), a.rows());

    } else {
      new (this) Base(
          ChainableStack::instance_->memalloc_.alloc_array<Scalar>(a.size()),
          a.rows(), a.cols());
    }
    Base::operator=(a);
    return *this;
  }
  /**
   * Forces hard copying matrices into an arena matrix
   * @tparam T Any type assignable to `Base`
   * @param x the values to write to `this`
   */
  template <typename T>
  void deep_copy(const T& x) {
    Base::operator=(x);
  }
};

}  // namespace math
}  // namespace stan

namespace Eigen {
namespace internal {

template <typename T>
struct traits<stan::math::arena_matrix<T>> {
  using base = traits<Eigen::Map<T>>;
  using XprKind = typename Eigen::internal::traits<std::decay_t<T>>::XprKind;
  enum {
    PlainObjectTypeInnerSize = base::PlainObjectTypeInnerSize,
    InnerStrideAtCompileTime = base::InnerStrideAtCompileTime,
    OuterStrideAtCompileTime = base::OuterStrideAtCompileTime,
    Alignment = base::Alignment,
    Flags = base::Flags
  };
};

}  // namespace internal
}  // namespace Eigen

#endif
