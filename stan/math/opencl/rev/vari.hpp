#ifndef STAN_MATH_OPENCL_REV_VARI_HPP
#define STAN_MATH_OPENCL_REV_VARI_HPP
#ifdef STAN_OPENCL

#include <stan/math/rev/core/vari.hpp>
#include <stan/math/opencl/is_matrix_cl.hpp>
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
   * Construct a dense matrix_cl variable implementation from a value. The
   * adjoint is initialized to zero.
   *
   * All constructed variables are added to the stack. Variables
   * should be constructed before variables on which they depend
   * to insure proper partial derivative propagation.  During
   * derivative propagation, the chain() method of each variable
   * will be called in the reverse order of construction.
   *
   * @tparam S A dense Eigen type that is convertible to `value_type`
   * @param x Value of the constructed variable.
   */
  template <typename S, require_convertible_t<S&, T>* = nullptr>
  explicit vari_value(S&& x)
      : chainable_alloc(),
        adj_(constant(0, x.rows(), x.cols())),
        val_(std::forward<S>(x)) {
    ChainableStack::instance_->var_stack_.push_back(this);
  }

 protected:
  /**
   * Construct a dense matrix_cl variable implementation from a value and
   * adjoint.
   *
   * All constructed variables are added to the stack. Variables
   * should be constructed before variables on which they depend
   * to insure proper partial derivative propagation.  During
   * derivative propagation, the chain() method of each variable
   * will be called in the reverse order of construction.
   *
   * @tparam R A dense Eigen type that is convertible to `T_val`
   * @tparam S A dense Eigen type that is convertible to `T_adj`
   * @param val Value of the constructed variable.
   * @param adj Adjoint of the constructed variable.
   */
  template <typename R, typename S, require_convertible_t<R&, T>* = nullptr,
            require_convertible_t<S&, T>* = nullptr>
  vari_value(R&& val, S&& adj)
      : chainable_alloc(),
        adj_(std::forward<S>(adj)),
        val_(std::forward<R>(val)) {
    ChainableStack::instance_->var_stack_.push_back(this);
  }

 public:
  /**
   * Construct an dense matrix_cl variable implementation from a value. The
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
   * @tparam S A dense Eigen type that is convertible to `value_type`
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
    const auto& val_block = stan::math::block(val_, row, col, rows, cols);
    const auto& adj_block = stan::math::block(adj_, row, col, rows, cols);
    return vari_value<std::decay_t<decltype(val_block)>>(val_block, adj_block);
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

  /**
   * Allocate memory from the underlying memory pool.  This memory is
   * is managed as a whole externally.
   *
   * Warning: Classes should not be allocated with this operator
   * if they have non-trivial destructors.
   *
   * @param nbytes Number of bytes to allocate.
   * @return Pointer to allocated bytes.
   */
  static inline void* operator new(size_t nbytes) {
    return ChainableStack::instance_->memalloc_.alloc(nbytes);
  }

  /**
   * Delete a pointer from the underlying memory pool.
   *
   * This no-op implementation enables a subclass to throw
   * exceptions in its constructor.  An exception thrown in the
   * constructor of a subclass will result in an error being
   * raised, which is in turn caught and calls delete().
   *
   * See the discussion of "plugging the memory leak" in:
   *   http://www.parashift.com/c++-faq/memory-pools.html
   */
  static inline void operator delete(void* /* ignore arg */) { /* no op */
  }

 private:
  template <typename>
  friend class var_value;
};

}  // namespace math
}  // namespace stan

#endif
#endif
