#ifndef STAN_MATH_OPENCL_REV_VARI_HPP
#define STAN_MATH_OPENCL_REV_VARI_HPP
#ifdef STAN_OPENCL

#include <stan/math/rev/core/vari.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/prim/reverse.hpp>

namespace stan {
namespace math {

template <typename T>
class vari_cl_base : public vari_base {
 public:
  /**
   * The value of this variable.
   */
  T val_;

  /**
   * The adjoint of this variable, which is the partial derivative
   * of this variable with respect to the root variable.
   */
  T adj_;

  /**
   * Construct a matrix_cl variable implementation from a value and
   * adjoint.
   *
   * This constructor does not add the matrix to any stack! It is intended only
   * for views into matrices already on stack.
   *
   * @param val Value of the constructed variable.
   * @param adj Adjoint of the constructed variable.
   */
  template <typename S, typename U,
            require_t<std::is_constructible<T, std::decay_t<S>>>* = nullptr,
            require_t<std::is_constructible<T, std::decay_t<U>>>* = nullptr>
  vari_cl_base(S&& val, U&& adj)
      : val_(std::forward<S>(val)), adj_(std::forward<U>(adj)) {}

  /**
   * Return a constant reference to the value of this vari.
   *
   * @return The value of this vari.
   */
  inline const auto& val() const { return val_; }
  inline auto& val_op() { return val_; }

  /**
   * Return a reference to the derivative of the root expression with
   * respect to this expression.  This method only works
   * after one of the `grad()` methods has been
   * called.
   *
   * @return Adjoint for this vari.
   */
  inline auto& adj() { return adj_; }
  inline auto& adj() const { return adj_; }
  inline auto& adj_op() { return adj_; }

  /**
   * Returns a view into a block of matrix.
   * @param row starting row of the block
   * @param col starting column of the block
   * @param rows number of rows in the block
   * @param cols number of columns in the block
   * @return block
   */
  auto block(int row, int col, int rows, int cols) {
    auto&& val_block = stan::math::block_zero_based(val_, row, col, rows, cols);
    auto&& adj_block = stan::math::block_zero_based(adj_, row, col, rows, cols);
    return vari_view<std::decay_t<decltype(val_block)>>(std::move(val_block),
                                                        std::move(adj_block));
  }

  /**
   * Returns a transposed view into the matrix.
   * @return transpose
   */
  auto transpose() {
    auto&& val_t = stan::math::transpose(val_);
    auto&& adj_t = stan::math::transpose(adj_);
    return vari_view<std::decay_t<decltype(val_t)>>(std::move(val_t),
                                                    std::move(adj_t));
  }

  /**
   * Returns column vector view into the row or column vector.
   * @return column vector view
   */
  auto as_column_vector_or_scalar() {
    auto&& val_t = stan::math::as_column_vector_or_scalar(val_);
    auto&& adj_t = stan::math::as_column_vector_or_scalar(adj_);
    return vari_view<std::decay_t<decltype(val_t)>>(std::move(val_t),
                                                    std::move(adj_t));
  }

  /**
   * Returns reverse view into the row or column vector.
   * @return reverse view
   */
  auto reverse() {
    auto&& val_t = stan::math::reverse(val_);
    auto&& adj_t = stan::math::reverse(adj_);
    return vari_view<std::decay_t<decltype(val_t)>>(std::move(val_t),
                                                    std::move(adj_t));
  }

  /**
   * Return indexed view into a matrix.
   *
   * Do not use with indices that reference any element of this more than once -
   * that cans cause data races in rev operations on the result!
   *
   * @param row_index kg expression used for row index
   * @param col_index kg expression used for column index
   */
  template <typename RowIndex, typename ColIndex>
  auto index(const RowIndex& row_index, const ColIndex& col_index) {
    RowIndex r1 = row_index;
    RowIndex r2 = row_index;
    ColIndex c1 = col_index;
    ColIndex c2 = col_index;
    auto&& val_t = stan::math::indexing(val_, std::move(r1), std::move(c1));
    auto&& adj_t = stan::math::indexing(adj_, std::move(r2), std::move(c2));
    return vari_view<std::decay_t<decltype(val_t)>>(std::move(val_t),
                                                    std::move(adj_t));
  }

  /**
   * Return the number of rows for this class's `val_` member
   */
  const Eigen::Index rows() const { return val_.rows(); }
  /**
   * Return the number of columns for this class's `val_` member
   */
  const Eigen::Index cols() const { return val_.cols(); }
  /**
   * Return the size of this class's `val_` member
   */
  const Eigen::Index size() const { return rows() * cols(); }

  virtual void chain() {}
};

template <typename T>
class vari_view<T, require_kernel_expression_lhs_t<T>>
    : public vari_cl_base<T> {
 public:
  /**
   * Rows at compile time
   */
  static constexpr int RowsAtCompileTime{-1};
  /**
   * Columns at compile time
   */
  static constexpr int ColsAtCompileTime{-1};

  using value_type = T;
  using vari_cl_base<T>::vari_cl_base;
  inline void set_zero_adjoint() final {}
};

/**
 * The variable implementation for `matrix_cl`.
 *
 * This class is complete (not abstract) and may be used for
 * constants.
 *
 * A variable implementation is constructed with a constant
 * value. It also stores the adjoint for storing the partial
 * derivative with respect to the root of the derivative tree.
 *
 */
template <typename T>
class vari_value<T, require_matrix_cl_t<T>> : public chainable_alloc,
                                              public vari_cl_base<T> {
 public:
  using value_type = T;

  /**
   * Rows at compile time
   */
  static constexpr int RowsAtCompileTime{-1};
  /**
   * Columns at compile time
   */
  static constexpr int ColsAtCompileTime{-1};

  /**
   * Construct a matrix_cl variable implementation from a value. The
   * adjoint is initialized to zero.
   *
   * All constructed variables are added to the stack. Variables
   * should be constructed before variables on which they depend
   * to insure proper partial derivative propagation.  During
   * derivative propagation, the chain() method of each variable
   * will be called in the reverse order of construction.
   *
   * @tparam S A `matrix_cl` or kernel generator expression type that is
   * convertible to `T`
   * @param x Value of the constructed variable.
   */
  template <typename S, require_convertible_t<S&, T>* = nullptr>
  explicit vari_value(S&& x)
      : chainable_alloc(),
        vari_cl_base<T>(std::forward<S>(x), constant(0, x.rows(), x.cols())) {
    ChainableStack::instance_->var_stack_.push_back(this);
  }

  /**
   * Construct a matrix_cl variable implementation from an Eigen value. The
   * adjoint is initialized to zero.
   *
   * All constructed variables are added to the stack. Variables
   * should be constructed before variables on which they depend
   * to insure proper partial derivative propagation.  During
   * derivative propagation, the chain() method of each variable
   * will be called in the reverse order of construction.
   *
   * @tparam S A dense Eigen value or expression that has same scalar type as T
   * @param x Value of the constructed variable.
   */
  template <typename S, require_eigen_t<S>* = nullptr,
            require_vt_same<T, S>* = nullptr>
  explicit vari_value(const S& x)
      : chainable_alloc(), vari_cl_base<T>(x, constant(0, x.rows(), x.cols())) {
    ChainableStack::instance_->var_nochain_stack_.push_back(this);
  }

  /**
   * Construct an matrix_cl variable implementation from a value. The
   *  adjoint is initialized to zero and if `stacked` is `false` this vari
   *  will be not be put on the var_stack. Instead it will only be put on
   *  a stack to keep track of whether the adjoint needs to be set to zero.
   *
   * All constructed variables are added to a stack.  Variables
   *  should be constructed before variables on which they depend
   *  to insure proper partial derivative propagation.  During
   *  derivative propagation, the chain() method of each variable
   *  will be called in the reverse order of construction.
   *
   * @tparam S A `matrix_cl` or kernel generator expression type that is
   * convertible to `T`
   * @param x Value of the constructed variable.
   * @param stacked If false will put this this vari on the nochain stack so
   * that its `chain()` method is not called.
   */
  template <typename S, require_convertible_t<S&, T>* = nullptr>
  vari_value(S&& x, bool stacked)
      : chainable_alloc(),
        vari_cl_base<T>(std::forward<S>(x), constant(0, x.rows(), x.cols())) {
    if (stacked) {
      ChainableStack::instance_->var_stack_.push_back(this);
    } else {
      ChainableStack::instance_->var_nochain_stack_.push_back(this);
    }
  }

  /**
   * Construct a dense Eigen variable implementation from a
   *  preconstructed values and adjoints.
   *
   * All constructed variables are not added to the stack. Variables
   * should be constructed before variables on which they depend
   * to insure proper partial derivative propagation.
   * @tparam S A dense Eigen type that is convertible to `value_type`
   * @tparam K A dense Eigen type that is convertible to `value_type`
   * @param val Matrix of values
   * @param adj Matrix of adjoints
   */
  template <typename S, typename K, require_convertible_t<T, S>* = nullptr,
            require_convertible_t<T, K>* = nullptr>
  explicit vari_value(S&& val, K&& adj)
      : chainable_alloc(),
        vari_cl_base<T>(std::forward<S>(val), std::forward<K>(adj)) {
    ChainableStack::instance_->var_nochain_stack_.push_back(this);
  }

  /**
   * Set the adjoint value of this variable to 0.  This is used to
   * reset adjoints before propagating derivatives again (for
   * example in a Jacobian calculation).
   */
  inline void set_zero_adjoint() final {
    this->adj_ = constant(0, this->rows(), this->cols());
  }

 private:
  template <typename, typename>
  friend class var_value;
};

}  // namespace math
}  // namespace stan

#endif
#endif
