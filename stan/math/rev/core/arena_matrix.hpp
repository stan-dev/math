#ifndef STAN_MATH_REV_CORE_ARENA_MATRIX_HPP
#define STAN_MATH_REV_CORE_ARENA_MATRIX_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core/chainable_alloc.hpp>
#include <stan/math/rev/core/chainable_object.hpp>
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

 private:
  /**
   * Move an rvalue value type'd Eigen matrix onto the `var_alloc_stack_`
   * @tparam T An rvalue Eigen matrix or vector
   * @param other whose underlying memory will be moved to an object sitting
   *  on Stan's `var_alloc_stack`.
   */
  template <typename T,
            std::enable_if_t<std::is_rvalue_reference<T&&>::value && is_plain_type<T>::value>* = nullptr>
  inline auto move_or_copy_map(T&& other) {
    using other_t = std::decay_t<T>;
    constexpr Eigen::Index OtherRowsAtCompileTime = other_t::RowsAtCompileTime;
    constexpr Eigen::Index OtherColsAtCompileTime = other_t::ColsAtCompileTime;
    auto&& chain_obj = make_chainable_ptr(std::forward<T>(other));
    return Base(
        const_cast<Scalar*>(chain_obj->data()),
        (RowsAtCompileTime == 1 && OtherColsAtCompileTime == 1)
                || (ColsAtCompileTime == 1 && OtherRowsAtCompileTime == 1)
            ? chain_obj->cols()
            : chain_obj->rows(),
        (RowsAtCompileTime == 1 && OtherColsAtCompileTime == 1)
                || (ColsAtCompileTime == 1 && OtherRowsAtCompileTime == 1)
            ? chain_obj->rows()
            : chain_obj->cols());
  }

  /**
   * Copy an lvalue type Eigen matrix into the arena memory
   * @tparam T An lvalue Eigen matrix or vector
   * @param other Eigen matrix or vector whose underlying memory is copied to Stan's arena allocator.
   */
  template <typename T,
            std::enable_if_t<!std::is_rvalue_reference<T&&>::value || !is_plain_type<T>::value>* = nullptr>
  inline auto move_or_copy_map(T&& other) {
    using other_t = std::decay_t<T>;
    constexpr Eigen::Index OtherRowsAtCompileTime = other_t::RowsAtCompileTime;
    constexpr Eigen::Index OtherColsAtCompileTime = other_t::ColsAtCompileTime;
    auto mapper = Base(
        ChainableStack::instance_->memalloc_.alloc_array<Scalar>(other.size()),
        (RowsAtCompileTime == 1 && OtherColsAtCompileTime == 1)
                || (ColsAtCompileTime == 1 && OtherRowsAtCompileTime == 1)
            ? other.cols()
            : other.rows(),
        (RowsAtCompileTime == 1 && OtherColsAtCompileTime == 1)
                || (ColsAtCompileTime == 1 && OtherRowsAtCompileTime == 1)
            ? other.rows()
            : other.cols());
    mapper = other;
    return mapper;
  }

 public:
  /**
   * Constructs `arena_matrix` from an eigen expression that is not derived from
   * a map.
   * @tparam T An Eigen type not inheriting from `Eigen::MapBase` that has
   * a `Scalar` type the is the same as `arena_matrix<MatrixType>::Scalar`.
   * @param other expression
   */
  template <typename T, require_eigen_t<T>* = nullptr,
            require_not_eigen_map_base_t<T>* = nullptr,
            require_same_t<Scalar, value_type_t<T>>* = nullptr>
  arena_matrix(T&& other)  // NOLINT
      : Base::Map(this->move_or_copy_map(std::forward<T>(other))) {}

  /**
   * Constructs `arena_matrix` from an eigen expression that is not derived from
   * a map.
   * @tparam T An Eigen type not inheriting from `Eigen::MapBase` that has
   * a `Scalar` type the is promotable to `arena_matrix<MatrixType>::Scalar`.
   * @param other expression
   */
  template <typename T, require_eigen_t<T>* = nullptr,
            require_not_eigen_map_base_t<T>* = nullptr,
            require_not_same_t<Scalar, value_type_t<T>>* = nullptr>
  arena_matrix(T&& other)  // NOLINT
      : Base::Map(
          ChainableStack::instance_->memalloc_.alloc_array<Scalar>(
              other.size()),
          (RowsAtCompileTime == 1 && std::decay_t<T>::ColsAtCompileTime == 1)
                  || (ColsAtCompileTime == 1
                      && std::decay_t<T>::RowsAtCompileTime == 1)
              ? other.cols()
              : other.rows(),
          (RowsAtCompileTime == 1 && std::decay_t<T>::ColsAtCompileTime == 1)
                  || (ColsAtCompileTime == 1
                      && std::decay_t<T>::RowsAtCompileTime == 1)
              ? other.rows()
              : other.cols()) {
    *this = other;
  }

  /**
   * Constructs `arena_matrix` from an expression. This makes an assumption that
   *  any other `Eigen::Map` contains memory that will exist for at least as
   *  long as the memory is needed for the reverse pass of reverse mode AD.
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
   * Constructs `arena_matrix` from an eigen expression that is not derived from
   * a map.
   * @tparam T An Eigen type not inheriting from `Eigen::MapBase` that has
   * a `Scalar` type the is the same as `arena_matrix<MatrixType>::Scalar`.
   * @param other expression
   */
  template <typename T, require_eigen_t<T>* = nullptr,
            require_not_eigen_map_base_t<T>* = nullptr,
            require_same_t<Scalar, value_type_t<T>>* = nullptr>
  arena_matrix& operator=(T&& other) {
    new (this) Base(this->move_or_copy_map(std::forward<T>(other)));
    return *this;
  }

  /**
   * Assignment operator for assigning an expression.
   * @param a expression to evaluate into this
   * @return `*this`
   */
  template <typename T, require_eigen_t<T>* = nullptr,
            require_not_eigen_map_base_t<T>* = nullptr,
            require_not_same_t<Scalar, value_type_t<T>>* = nullptr>
  arena_matrix& operator=(T&& a) {
    // do we need to transpose?
    if ((RowsAtCompileTime == 1 && std::decay_t<T>::ColsAtCompileTime == 1)
        || (ColsAtCompileTime == 1
            && std::decay_t<T>::RowsAtCompileTime == 1)) {
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
