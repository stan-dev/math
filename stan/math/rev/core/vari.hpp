#ifndef STAN_MATH_REV_CORE_VARI_HPP
#define STAN_MATH_REV_CORE_VARI_HPP

#include <stan/math/rev/core/chainable_alloc.hpp>
#include <stan/math/rev/core/chainablestack.hpp>
#include <stan/math/prim/meta.hpp>
#include <ostream>
#include <type_traits>

namespace stan {
namespace math {

// forward declaration of var
template <typename T, typename>
class var_value;

/**
 * The variable implementation base class.
 *
 * This class is complete (not abstract) and may be used for
 * constants.
 *
 * A variable implementation is constructed with a constant
 * value. It also stores the adjoint for storing the partial
 * derivative with respect to the root of the derivative tree.
 *
 * The chain() method applies the chain rule. Concrete extensions
 * of this class will represent base variables or the result
 * of operations such as addition or subtraction. These extended
 * classes will store operand variables and propagate derivative
 * information via an implementation of chain().
 */
template <typename T, typename>
class vari_value;

template <typename T>
class vari_value<T, std::enable_if_t<std::is_floating_point<T>::value>> {
 private:
  template <typename, typename>
  friend class var_value;

 public:
  using Scalar = T;
  using value_type = Scalar;
  /**
   * The value of this variable.
   */
  const Scalar val_;

  /**
   * The adjoint of this variable, which is the partial derivative
   * of this variable with respect to the root variable.
   */
  Scalar adj_;

  /**
   * Construct a variable implementation from a value.  The
   * adjoint is initialized to zero.
   *
   * All constructed variables are added to the stack.  Variables
   * should be constructed before variables on which they depend
   * to insure proper partial derivative propagation.  During
   * derivative propagation, the chain() method of each variable
   * will be called in the reverse order of construction.
   *
   * @tparam S an Arithmetic type.
   * @param x Value of the constructed variable.
   */
  template <typename S,
            std::enable_if_t<std::is_convertible<S&, Scalar>::value>* = nullptr>
  vari_value(S x) noexcept : val_(x), adj_(0.0) {  // NOLINT
    ChainableStack::instance_->var_stack_.emplace_back(this);
  }

  /**
   * Construct a variable implementation from a value.  The
   *  adjoint is initialized to zero and if `stacked` is `false` this vari
   *  will be moved to the nochain stack s.t. it's chain method will not be
   *  called when calling `grad()`.
   *
   * All constructed variables are added to the stack.  Variables
   *  should be constructed before variables on which they depend
   *  to insure proper partial derivative propagation.  During
   *  derivative propagation, the chain() method of each variable
   *  will be called in the reverse order of construction.
   *
   * @tparam S an Arithmetic type.
   * @param x Value of the constructed variable.
   * @param stacked If false will put this this vari on the nochain stack so
   * that its `chain()` method is not called.
   */
  template <typename S,
            std::enable_if_t<std::is_convertible<S&, Scalar>::value>* = nullptr>
  vari_value(S x, bool stacked) noexcept : val_(x), adj_(0.0) {
    if (stacked) {
      ChainableStack::instance_->var_stack_.emplace_back(this);
    } else {
      ChainableStack::instance_->var_nochain_stack_.emplace_back(this);
    }
  }

  /**
   * Constructor from vari_value
   * @tparam S An arithmetic type
   * @param x A vari_value
   */
  vari_value(const vari_value<T>& x) noexcept : val_(x.val_), adj_(x.adj_) {
    ChainableStack::instance_->var_stack_.emplace_back(this);
  }
  /**
   * Constructor from vari_value
   * @tparam S An arithmetic type
   * @param x A vari_value
   */
  vari_value(vari_value<T>&& x) noexcept : val_(x.val_), adj_(x.adj_) {
    ChainableStack::instance_->var_stack_.emplace_back(this);
  }

  /**
   * Initialize the adjoint for this (dependent) variable to 1.
   * This operation is applied to the dependent variable before
   * propagating derivatives, setting the derivative of the
   * result with respect to itself to be 1.
   */
  inline void init_dependent() noexcept { adj_ = 1.0; }

  /**
   * Set the adjoint value of this variable to 0.  This is used to
   * reset adjoints before propagating derivatives again (for
   * example in a Jacobian calculation).
   */
  inline void set_zero_adjoint() noexcept { adj_ = 0.0; }

  virtual void chain() {}

  /**
   * Insertion operator for vari. Prints the current value and
   * the adjoint value.
   *
   * @param os [in, out] ostream to modify
   * @param v [in] vari object to print.
   *
   * @return The modified ostream.
   */
  friend std::ostream& operator<<(std::ostream& os, const vari_value<T>* v) {
    return os << v->val_ << ":" << v->adj_;
  }

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
  static inline void* operator new(size_t nbytes) noexcept {
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
  static inline void operator delete(
      void* /* ignore arg */) noexcept { /* no op */
  }
};

// For backwards compatability the default is double
using vari = vari_value<double>;



template <typename T>
class vari_value<T, std::enable_if_t<is_eigen<T>::value>> {
 private:
  template <typename, typename>
  friend class var_value;
  using eigen_scalar = value_type_t<T>;
  value_type_t<T>* val_mem_;
  value_type_t<T>* adj_mem_;

 public:
  using Scalar = T;
  using value_type = T;
  using PlainObject = typename T::PlainObject;
  using eigen_map = Eigen::Map<PlainObject>;
  const Eigen::Index rows_;
  const Eigen::Index cols_;
  const Eigen::Index size_;
  const Eigen::Index RowsAtCompileTime = T::RowsAtCompileTime;
  const Eigen::Index ColsAtCompileTime = T::ColsAtCompileTime;
  /**
   * The value of this variable.
   */
  eigen_map val_;

  /**
   * The adjoint of this variable, which is the partial derivative
   * of this variable with respect to the root variable.
   */
  eigen_map adj_;

  /**
   * Construct a variable implementation from a value.  The
   * adjoint is initialized to zero.
   *
   * All constructed variables are added to the stack.  Variables
   * should be constructed before variables on which they depend
   * to insure proper partial derivative propagation.  During
   * derivative propagation, the chain() method of each variable
   * will be called in the reverse order of construction.
   *
   * @param x Value of the constructed variable.
   */
  explicit vari_value(const T& x)
      : val_mem_(ChainableStack::instance_->memalloc_.alloc_array<eigen_scalar>(
            x.size())),
        adj_mem_(ChainableStack::instance_->memalloc_.alloc_array<eigen_scalar>(
            x.size())),
        rows_(x.rows()),
        cols_(x.cols()),
        size_(x.size()),
        val_(val_mem_, x.rows(), x.cols()),
        adj_(adj_mem_, x.rows(), x.cols()) {
    val_ = x;
    adj_.setZero();
    ChainableStack::instance_->var_stack_.push_back(this);
  }

  vari_value(const T& x, bool stacked)
      : val_mem_(ChainableStack::instance_->memalloc_.alloc_array<eigen_scalar>(
            x.size())),
        adj_mem_(ChainableStack::instance_->memalloc_.alloc_array<eigen_scalar>(
            x.size())),
        rows_(x.rows()),
        cols_(x.cols()),
        size_(x.size()),
        val_(val_mem_, x.rows(), x.cols()),
        adj_(adj_mem_, x.rows(), x.cols()) {
    val_ = x;
    adj_.setZero();
    if (stacked) {
      ChainableStack::instance_->var_stack_.push_back(this);
    } else {
      ChainableStack::instance_->var_nochain_stack_.push_back(this);
    }
  }

  const Eigen::Index rows() const { return rows_; }
  const Eigen::Index cols() const { return cols_; }
  const Eigen::Index size() const { return size_; }
  virtual void chain() {}
  /**
   * Initialize the adjoint for this (dependent) variable to 1.
   * This operation is applied to the dependent variable before
   * propagating derivatives, setting the derivative of the
   * result with respect to itself to be 1.
   */
  void init_dependent() { adj_.setOnes(); }

  /**
   * Set the adjoint value of this variable to 0.  This is used to
   * reset adjoints before propagating derivatives again (for
   * example in a Jacobian calculation).
   */
  void set_zero_adjoint() { adj_.setZero(); }

  /**
   * Insertion operator for vari. Prints the current value and
   * the adjoint value.
   *
   * @param os [in, out] ostream to modify
   * @param v [in] vari object to print.
   *
   * @return The modified ostream.
   */
  friend std::ostream& operator<<(std::ostream& os, const vari_value<T>* v) {
    return os << v->val_ << ":" << v->adj_;
  }

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
};

struct vari_printer {
  std::ostream& o_;
  int i_{0};
  vari_printer(std::ostream& o, int i) : o_(o), i_(i) {}
  template <typename T>
  inline void operator()(T*& x) const {
    o_ << i_ << "  " << x << "  " << x->val_ << " : " << x->adj_ << std::endl;
  }
};
}  // namespace math
}  // namespace stan
#endif
