#ifndef STAN_MATH_OPENCL_REV_VARI_HPP
#define STAN_MATH_OPENCL_REV_VARI_HPP
#ifdef STAN_OPENCL

#include <stan/math/rev/core/vari.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/opencl/kernel_generator.hpp>

namespace stan {
namespace math {
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
class vari_value<T, require_kernel_expression_lhs_t<T>>
    : public vari_base, public chainable_alloc {
 public:
  using value_type = T;
  /**
   * The adjoint of this variable, which is the partial derivative
   * of this variable with respect to the root variable.
   */
  T adj_;

  /**
   * The value of this variable.
   */
  T val_;

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
        adj_(constant(0, x.rows(), x.cols())),
        val_(std::forward<S>(x)) {
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
      : chainable_alloc(), adj_(constant(0, x.rows(), x.cols())), val_(x) {
    ChainableStack::instance_->var_stack_.push_back(this);
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
        adj_(constant(0, x.rows(), x.cols())),
        val_(std::forward<S>(x)) {
    if (stacked) {
      ChainableStack::instance_->var_stack_.push_back(this);
    } else {
      ChainableStack::instance_->var_nochain_stack_.push_back(this);
    }
  }

  /**
   * Returns a view into a block of matrix.
   * @param row starting row of the block
   * @param col starting column of the block
   * @param rows number of rows in the block
   * @param cols number of columns in the block
   * @return block
   */
  auto block(int row, int col, int rows, int cols) {
    auto&& val_block = stan::math::block(val_, row, col, rows, cols);
    auto&& adj_block = stan::math::block(adj_, row, col, rows, cols);
    return vari_value<std::decay_t<decltype(val_block)>>(std::move(val_block),
                                                         std::move(adj_block));
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
  const Eigen::Index size() const { return val_.size(); }

  virtual void chain() {}

  /**
   * Set the adjoint value of this variable to 0.  This is used to
   * reset adjoints before propagating derivatives again (for
   * example in a Jacobian calculation).
   */
  inline void set_zero_adjoint() final { adj_ = constant(0, rows(), cols()); }

 protected:
  // to allow access to this constructor from instantinations with different
  // template parameters
  template <typename, typename>
  friend class vari_value;

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
  vari_value(T&& val, T&& adj)
      : chainable_alloc(),
        adj_(std::forward<T>(adj)),
        val_(std::forward<T>(val)) {}

 private:
  template <typename, typename>
  friend class var_value;
};

}  // namespace math
}  // namespace stan

#endif
#endif
