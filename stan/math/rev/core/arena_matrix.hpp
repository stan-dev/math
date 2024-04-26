#ifndef STAN_MATH_REV_CORE_ARENA_MATRIX_HPP
#define STAN_MATH_REV_CORE_ARENA_MATRIX_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/core/chainable_alloc.hpp>
#include <stan/math/rev/core/chainablestack.hpp>
#include <stan/math/rev/core/chainable_object.hpp>
#include <stan/math/rev/core/var_value_fwd_declare.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
namespace stan {
namespace math {

template <typename MatrixType>
class arena_matrix<MatrixType, require_eigen_dense_base_t<MatrixType>>
    : public Eigen::Map<MatrixType> {
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
            require_plain_type_t<T>* = nullptr,
            require_same_t<T, MatrixType>* = nullptr>
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
            require_plain_type_t<T>* = nullptr,
            require_same_t<T, MatrixType>* = nullptr>
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

template <typename MatrixType>
class arena_matrix<MatrixType, require_eigen_sparse_base_t<MatrixType>>
    : public Eigen::Map<MatrixType> {
 public:
  using Scalar = value_type_t<MatrixType>;
  using Base = Eigen::Map<MatrixType>;
  using PlainObject = std::decay_t<MatrixType>;
  using StorageIndex = typename PlainObject::StorageIndex;

  /**
   * Default constructor.
   */
  arena_matrix() : Base::Map(0, 0, 0, nullptr, nullptr, nullptr) {}
  /**
   * Constructor for CDC formatted data. Assumes compressed and does not own
   * memory.
   *
   * Constructs a read-write Map to a sparse matrix of size rows x cols,
   * containing nnz non-zero coefficients, stored as a sparse format as defined
   * by the pointers outerIndexPtr, innerIndexPtr, and valuePtr. If the optional
   * parameter innerNonZerosPtr is the null pointer, then a standard compressed
   * format is assumed.
   *
   * @param rows Number of rows
   * @param cols Number of columns
   * @param nnz Number of nonzero items in matrix
   * @param outerIndexPtr A pointer to the array of outer indices.
   * @param innerIndexPtr A pointer to the array of inner indices.
   * @param valuePtr A pointer to the array of values.
   * @param innerNonZerosPtr A pointer to the array of the number of non zeros
   * of the inner vectors.
   *
   */
  arena_matrix(Eigen::Index rows, Eigen::Index cols, Eigen::Index nnz,
               StorageIndex* outerIndexPtr, StorageIndex* innerIndexPtr,
               Scalar* valuePtr, StorageIndex* innerNonZerosPtr = nullptr)
      : Base::Map(rows, cols, nnz, outerIndexPtr, innerIndexPtr, valuePtr,
                  innerNonZerosPtr) {}

 private:
  template <typename T, typename Integral>
  inline T* copy_vector(const T* ptr, Integral size) {
    if (size == 0)
      return nullptr;
    T* ret = ChainableStack::instance_->memalloc_.alloc_array<T>(size);
    std::copy_n(ptr, size, ret);
    return ret;
  }

 public:
  /**
   * Constructs `arena_matrix` from an `Eigen::SparseMatrix`.
   * @param other Eigen Sparse Matrix class
   */
  template <typename T, require_same_t<T, PlainObject>* = nullptr>
  arena_matrix(T&& other)  // NOLINT
      : Base::Map(
          other.rows(), other.cols(), other.nonZeros(),
          copy_vector(other.outerIndexPtr(), other.outerSize() + 1),
          copy_vector(other.innerIndexPtr(), other.nonZeros()),
          copy_vector(other.valuePtr(), other.nonZeros()),
          copy_vector(
              other.innerNonZeroPtr(),
              other.innerNonZeroPtr() == nullptr ? 0 : other.innerSize())) {}

  /**
   * Constructs `arena_matrix` from an Eigen expression
   * @param other An expression
   */
  template <typename S, require_convertible_t<S&, PlainObject>* = nullptr,
            require_not_same_t<S, PlainObject>* = nullptr>
  arena_matrix(S&& other)  // NOLINT
      : arena_matrix(PlainObject(std::forward<S>(other))) {}

  /**
   * Constructs `arena_matrix` from an expression. This makes an assumption that
   * any other `Eigen::Map` also contains memory allocated in the arena.
   * @param other expression
   */
  arena_matrix(const Base& other)  // NOLINT
      : Base(other) {}

  /**
   * Const copy constructor. No actual copy is performed
   * @note Since the memory for the arena matrix sits in Stan's memory arena all
   * copies/moves of arena matrices are shallow moves
   * @param other matrix to copy from
   */
  arena_matrix(const arena_matrix<MatrixType>& other)
      : Base::Map(other.rows(), other.cols(), other.nonZeros(),
                  const_cast<StorageIndex*>(other.outerIndexPtr()),
                  const_cast<StorageIndex*>(other.innerIndexPtr()),
                  const_cast<Scalar*>(other.valuePtr()),
                  const_cast<StorageIndex*>(other.innerNonZeroPtr())) {}
  /**
   * Move constructor.
   * @note Since the memory for the arena matrix sits in Stan's memory arena all
   * copies/moves of arena matrices are shallow moves
   * @param other matrix to copy from
   */
  arena_matrix(arena_matrix<MatrixType>&& other)
      : Base::Map(other.rows(), other.cols(), other.nonZeros(),
                  const_cast<StorageIndex*>(other.outerIndexPtr()),
                  const_cast<StorageIndex*>(other.innerIndexPtr()),
                  const_cast<Scalar*>(other.valuePtr()),
                  const_cast<StorageIndex*>(other.innerNonZeroPtr())) {}
  /**
   * Copy constructor. No actual copy is performed
   * @note Since the memory for the arena matrix sits in Stan's memory arena all
   * copies/moves of arena matrices are shallow moves
   * @param other matrix to copy from
   */
  arena_matrix(arena_matrix<MatrixType>& other)
      : Base::Map(other.rows(), other.cols(), other.nonZeros(),
                  const_cast<StorageIndex*>(other.outerIndexPtr()),
                  const_cast<StorageIndex*>(other.innerIndexPtr()),
                  const_cast<Scalar*>(other.valuePtr()),
                  const_cast<StorageIndex*>(other.innerNonZeroPtr())) {}

  // without this using, compiler prefers combination of implicit construction
  // and copy assignment to the inherited operator when assigned an expression
  using Base::operator=;

  /**
   * Copy assignment operator.
   * @tparam ArenaMatrix An `arena_matrix` type
   * @param other matrix to copy from
   * @return `*this`
   */
  template <typename ArenaMatrix,
            require_same_t<std::decay_t<ArenaMatrix>,
                           arena_matrix<MatrixType>>* = nullptr>
  arena_matrix& operator=(ArenaMatrix&& other) {
    // placement new changes what data map points to - there is no allocation
    new (this) Base(other.rows(), other.cols(), other.nonZeros(),
                    const_cast<StorageIndex*>(other.outerIndexPtr()),
                    const_cast<StorageIndex*>(other.innerIndexPtr()),
                    const_cast<Scalar*>(other.valuePtr()),
                    const_cast<StorageIndex*>(other.innerNonZeroPtr()));
    return *this;
  }

  /**
   * Assignment operator for assigning an expression.
   * @tparam Expr An expression that an `arena_matrix<MatrixType>` can be
   * constructed from
   * @param expr expression to evaluate into this
   * @return `*this`
   */
  template <typename Expr,
            require_not_same_t<Expr, arena_matrix<MatrixType>>* = nullptr>
  arena_matrix& operator=(Expr&& expr) {
    new (this) arena_matrix(std::forward<Expr>(expr));
    return *this;
  }

 private:
  /**
   * inplace operations functor for `Sparse (.*)= Sparse`.
   * @note This assumes that each matrix is of the same size and sparsity
   * pattern.
   * @tparam F A type with a valid `operator()(Scalar& x, const Scalar& y)`
   * method
   * @tparam Expr Type derived from `Eigen::SparseMatrixBase`
   * @param f Functor that performs the inplace operation
   * @param other The right hand side of the inplace operation
   */
  template <typename F, typename Expr,
            require_convertible_t<Expr&, MatrixType>* = nullptr,
            require_same_t<Expr, arena_matrix<MatrixType>>* = nullptr>
  inline void inplace_ops_impl(F&& f, Expr&& other) {
    auto&& x = to_ref(other);
    auto* val_ptr = (*this).valuePtr();
    auto* x_val_ptr = x.valuePtr();
    const auto non_zeros = (*this).nonZeros();
    for (Eigen::Index i = 0; i < non_zeros; ++i) {
      f(val_ptr[i], x_val_ptr[i]);
    }
  }

  /**
   * inplace operations functor for `Sparse (.*)= Sparse`.
   * @note This assumes that each matrix is of the same size and sparsity
   * pattern.
   * @tparam F A type with a valid `operator()(Scalar& x, const Scalar& y)`
   * method
   * @tparam Expr Type derived from `Eigen::SparseMatrixBase`
   * @param f Functor that performs the inplace operation
   * @param other The right hand side of the inplace operation
   */
  template <typename F, typename Expr,
            require_convertible_t<Expr&, MatrixType>* = nullptr,
            require_not_same_t<Expr, arena_matrix<MatrixType>>* = nullptr>
  inline void inplace_ops_impl(F&& f, Expr&& other) {
    auto&& x = to_ref(other);
    for (int k = 0; k < (*this).outerSize(); ++k) {
      typename Base::InnerIterator it(*this, k);
      typename std::decay_t<Expr>::InnerIterator iz(x, k);
      for (; static_cast<bool>(it) && static_cast<bool>(iz); ++it, ++iz) {
        f(it.valueRef(), iz.value());
      }
    }
  }

  /**
   * inplace operations functor for `Sparse (.*)= Dense`.
   * @note This assumes the user intends to perform the inplace operation for
   * the nonzero parts of `this`
   * @tparam F A type with a valid `operator()(Scalar& x, const Scalar& y)`
   * method
   * @tparam Expr Type derived from `Eigen::DenseBase`
   * @param f Functor that performs the inplace operation
   * @param other The right hand side of the inplace operation
   */
  template <typename F, typename Expr,
            require_not_convertible_t<Expr&, MatrixType>* = nullptr,
            require_eigen_dense_base_t<Expr>* = nullptr>
  inline void inplace_ops_impl(F&& f, Expr&& other) {
    auto&& x = to_ref(other);
    for (int k = 0; k < (*this).outerSize(); ++k) {
      typename Base::InnerIterator it(*this, k);
      for (; static_cast<bool>(it); ++it) {
        f(it.valueRef(), x(it.row(), it.col()));
      }
    }
  }

  /**
   * inplace operations functor for `Sparse (.*)= Scalar`.
   * @note This assumes the user intends to perform the inplace operation for
   * the nonzero parts of `this`
   * @tparam F A type with a valid `operator()(Scalar& x, const Scalar& y)`
   * method
   * @tparam T A scalar type
   * @param f Functor that performs the inplace operation
   * @param other The right hand side of the inplace operation
   */
  template <typename F, typename T,
            require_convertible_t<T&, Scalar>* = nullptr>
  inline void inplace_ops_impl(F&& f, T&& other) {
    auto* val_ptr = (*this).valuePtr();
    const auto non_zeros = (*this).nonZeros();
    for (Eigen::Index i = 0; i < non_zeros; ++i) {
      f(val_ptr[i], other);
    }
  }

 public:
  /**
   * Inplace operator `+=`
   * @note Caution! Inplace operators assume that either
   *  1. The right hand side sparse matrix has the same sparsity pattern
   *  2. You only intend to add a scalar or dense matrix coefficients to the
   * nonzero values of `this`. This is intended to be used within the reverse
   * pass for accumulation of the adjoint and is built as such. Any other use
   * case should be be sure that the above assumptions are satisfied.
   * @tparam T A type derived from `Eigen::SparseMatrixBase` or
   * `Eigen::DenseMatrixBase` or a `Scalar`
   * @param other value to be added inplace to the matrix.
   */
  template <typename T>
  inline arena_matrix& operator+=(T&& other) {
    inplace_ops_impl(
        [](auto&& x, auto&& y) {
          x += y;
          return;
        },
        std::forward<T>(other));
    return *this;
  }

  /**
   * Inplace operator `+=`
   * @note Caution! Inplace operators assume that either
   *  1. The right hand side sparse matrix has the same sparsity pattern
   *  2. You only intend to add a scalar or dense matrix coefficients to the
   * nonzero values of `this`. This is intended to be used within the reverse
   * pass for accumulation of the adjoint and is built as such. Any other use
   * case should be be sure that the above assumptions are satisfied.
   * @tparam T A type derived from `Eigen::SparseMatrixBase` or
   * `Eigen::DenseMatrixBase` or a `Scalar`
   * @param other value to be added inplace to the matrix.
   */
  template <typename T>
  inline arena_matrix& operator-=(T&& other) {
    inplace_ops_impl(
        [](auto&& x, auto&& y) {
          x -= y;
          return;
        },
        std::forward<T>(other));
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
