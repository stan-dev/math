#ifndef STAN_MATH_REV_CORE_ARENA_MATRIX_HPP
#define STAN_MATH_REV_CORE_ARENA_MATRIX_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/core/chainable_alloc.hpp>
#include <stan/math/rev/core/chainablestack.hpp>
#include <stan/math/rev/core/chainable_object.hpp>
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

 private:
  template <typename T>
  constexpr auto get_rows(const T& x) {
    return (RowsAtCompileTime == 1 && T::ColsAtCompileTime == 1)
                   || (ColsAtCompileTime == 1 && T::RowsAtCompileTime == 1)
               ? x.cols()
               : x.rows();
  }
  template <typename T>
  constexpr auto get_cols(const T& x) {
    return (RowsAtCompileTime == 1 && T::ColsAtCompileTime == 1)
                   || (ColsAtCompileTime == 1 && T::RowsAtCompileTime == 1)
               ? x.rows()
               : x.cols();
  }

 public:
  /**
   * Constructs `arena_matrix` from an expression
   * @param other expression
   */
  template <typename T, require_eigen_t<T>* = nullptr>
  arena_matrix(const T& other)  // NOLINT
      : Base::Map(ChainableStack::instance_->memalloc_.alloc_array<Scalar>(
                      other.size()),
                  get_rows(other), get_cols(other)) {
    *this = other;
  }
  /**
   * Overwrite the current arena_matrix with new memory and assign a matrix to
   * it
   * @tparam T An eigen type inheriting from `Eigen::EigenBase`
   * @param other A matrix that will be copied over to the arena allocator
   */
  template <typename T, require_eigen_t<T>* = nullptr>
  arena_matrix& operator=(const T& other) {
    new (this) Base(
        ChainableStack::instance_->memalloc_.alloc_array<Scalar>(other.size()),
        get_rows(other), get_cols(other));
    Base::operator=(other);
    return *this;
  }

  /**
   * Constructs `arena_matrix` from an rvalue expression that is a `plain_type`,
   *  then movies it to the object stack.
   * @tparam T A type that inherits from Eigen::DenseBase that is not an
   * `arena_matrix`.
   * @param other expression
   * @note When T is both an rvalue and a plain type, the expression is moved to
   * the object stack.
   */
  template <typename T, require_eigen_t<T>* = nullptr,
            require_not_arena_matrix_t<T>* = nullptr,
            require_t<std::is_rvalue_reference<T&&>>* = nullptr,
            require_plain_type_t<T>* = nullptr>
  arena_matrix(T&& other)  // NOLINT
      : Base::Map([](auto&& x) {
          using base_map_t =
              typename stan::math::arena_matrix<MatrixType>::Base;
          auto other_ptr = make_chainable_ptr(std::move(x));
          // other has it's rows and cols swapped already if it needed that
          return base_map_t(&(other_ptr->coeffRef(0)), other_ptr->rows(),
                            other_ptr->cols());
        }(std::move(other))) {}

  /**
   * Assignment operator for assigning an expression.
   * This is for rvalue plain type objects that can be moved over to the object
   * stack instead of allocating new memory.
   * @param other expression to evaluate into this
   * @return `*this`
   */
  template <typename T, require_eigen_t<T>* = nullptr,
            require_not_arena_matrix_t<T>* = nullptr,
            require_t<std::is_rvalue_reference<T&&>>* = nullptr,
            require_plain_type_t<T>* = nullptr>
  arena_matrix& operator=(T&& other) {
    auto other_ptr = make_chainable_ptr(std::move(other));
    new (this)
        Base(&(other_ptr->coeffRef(0)), other_ptr->rows(), other_ptr->cols());
    return *this;
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
