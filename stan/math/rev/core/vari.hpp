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
template <typename T>
class var_value;
/**
 * Abstract base class that all `vari_value` and it's derived classes inherit.
 *
 * The chain() method applies the chain rule. Concrete extensions
 * of this class will represent base variables or the result
 * of operations such as addition or subtraction. These extended
 * classes will store operand variables and propagate derivative
 * information via an implementation of chain().
 */
class vari_base {
 public:
  /**
   * Apply the chain rule to this variable based on the variables
   * on which it depends.
   */
  virtual void chain() = 0;
  virtual void set_zero_adjoint() = 0;
  virtual ~vari_base() noexcept {}
};

template <typename T, typename = void>
class vari_value;
/**
 * The variable implementation for floating point types.
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
class vari_value<T, require_floating_point_t<T>> : public vari_base {
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
   * @tparam S a floating point type.
   * @param x Value of the constructed variable.
   */
  template <typename S, require_convertible_t<S&, Scalar>* = nullptr>
  vari_value(S x) noexcept : val_(x), adj_(0.0) {  // NOLINT
    ChainableStack::instance_->var_stack_.emplace_back(this);
  }

  /**
   * Construct a variable implementation from a value.  The
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
   * @tparam S n floating point type.
   * @param x Value of the constructed variable.
   * @param stacked If false will put this this vari on the nochain stack so
   * that its `chain()` method is not called.
   */
  template <typename S, require_convertible_t<S&, Scalar>* = nullptr>
  vari_value(S x, bool stacked) noexcept : val_(x), adj_(0.0) {
    if (stacked) {
      ChainableStack::instance_->var_stack_.emplace_back(this);
    } else {
      ChainableStack::instance_->var_nochain_stack_.emplace_back(this);
    }
  }

  ~vari_value() = default;

  inline void chain() {}

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
  inline void set_zero_adjoint() noexcept final { adj_ = 0.0; }

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

 private:
  template <typename>
  friend class var_value;
};

// For backwards compatability the default is double
using vari = vari_value<double>;

/**
 * The variable implementation for Eigen dense matrix types.
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
class vari_value<T, require_eigen_dense_base_t<T>> : public vari_base {
 public:
  static_assert(
      is_plain_type<T>::value,
      "The template for this var is an"
      " expression but a var_value's inner type must be assignable such as"
      " a double, Eigen::Matrix, or Eigen::Array");

  /**
   * `PlainObject` represents a user constructible type such as Matrix or Array
   */
  using PlainObject = std::decay_t<plain_type_t<T>>;
  using Scalar = PlainObject;  // The underlying type for this class
  using value_type = Scalar;   // The underlying type for this class
  using eigen_scalar = value_type_t<PlainObject>;  // A floating point type
  /**
   * Maps for adj_ and val_
   */
  using eigen_map = Eigen::Map<PlainObject, Eigen::Aligned8>;
  using vari_type = vari_value<T, require_eigen_dense_base_t<T>>;
  eigen_scalar* val_mem_;  // Pointer to memory allocated on the stack for val_
  eigen_scalar* adj_mem_;  // Pointer to memory allocated on the stack for adj_
  /**
   * Number of rows known at compile time
   */
  static constexpr Eigen::Index RowsAtCompileTime = Scalar::RowsAtCompileTime;
  /**
   * Number of columns known at compile time
   */
  static constexpr Eigen::Index ColsAtCompileTime = Scalar::ColsAtCompileTime;

  /**
   * The value of this variable.
   */
  const eigen_map val_;

  /**
   * The adjoint of this variable, which is the partial derivative
   * of this variable with respect to the root variable.
   */
  eigen_map adj_;

  /**
   * Construct a dense Eigen variable implementation from a value. The
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
  template <typename S, require_convertible_t<S&, value_type>* = nullptr>
  explicit vari_value(S&& x)
      : val_mem_(ChainableStack::instance_->memalloc_.alloc_array<eigen_scalar>(
            x.size())),
        adj_mem_(ChainableStack::instance_->memalloc_.alloc_array<eigen_scalar>(
            x.size())),
        val_(eigen_map(val_mem_, x.rows(), x.cols()) = x),
        adj_(eigen_map(adj_mem_, x.rows(), x.cols()).setZero()) {
    ChainableStack::instance_->var_stack_.push_back(this);
  }

  /**
   * Construct an dense Eigen variable implementation from a value. The
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
  template <typename S, require_convertible_t<S&, value_type>* = nullptr>
  vari_value(S&& x, bool stacked)
      : val_mem_(ChainableStack::instance_->memalloc_.alloc_array<eigen_scalar>(
            x.size())),
        adj_mem_(ChainableStack::instance_->memalloc_.alloc_array<eigen_scalar>(
            x.size())),
        val_(eigen_map(val_mem_, x.rows(), x.cols()) = x),
        adj_(eigen_map(adj_mem_, x.rows(), x.cols()).setZero()) {
    if (stacked) {
      ChainableStack::instance_->var_stack_.push_back(this);
    } else {
      ChainableStack::instance_->var_nochain_stack_.push_back(this);
    }
  }

  /**
   * Return the number of rows for this class's `val_` member
   */
  const Eigen::Index rows() const { return val_.rows(); }
  /**
   * Return the number of columns for this class's `val_` member
   */
  const Eigen::Index cols() const { return val_.rows(); }
  /**
   * Return the size of this class's `val_` member
   */
  const Eigen::Index size() const { return val_.size(); }

  virtual void chain() {}
  /**
   * Initialize the adjoint for this (dependent) variable to 1.
   * This operation is applied to the dependent variable before
   * propagating derivatives, setting the derivative of the
   * result with respect to itself to be 1.
   */
  inline void init_dependent() { adj_.setOnes(); }

  /**
   * Set the adjoint value of this variable to 0.  This is used to
   * reset adjoints before propagating derivatives again (for
   * example in a Jacobian calculation).
   */
  inline void set_zero_adjoint() final { adj_.setZero(); }

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
    return os << "val: \n" << v->val_ << " \nadj: \n" << v->adj_;
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

 private:
  template <typename>
  friend class var_value;
  template <typename, typename>
  friend class vari_value;
};

/**
 * The variable implementation for Eigen sparse matrix types.
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
class vari_value<T, std::enable_if_t<is_eigen_sparse_base<T>::value>>
    : public vari_base, chainable_alloc {
 public:
  using PlainObject
      = std::decay_t<plain_type_t<T>>;  // Base type of Eigen class
  using Scalar = PlainObject;           // vari's adj_ and val_ member type
  using value_type = Scalar;            // vari's adj_ and val_ member type
  /**
   * Rows at compile time
   */
  static constexpr Eigen::Index RowsAtCompileTime = T::RowsAtCompileTime;
  /**
   * Columns at compile time
   */
  static constexpr Eigen::Index ColsAtCompileTime = T::ColsAtCompileTime;
  /**
   * The value of this variable.
   */
  const PlainObject val_;

  /**
   * The adjoint of this variable, which is the partial derivative
   * of this variable with respect to the root variable.
   */
  PlainObject adj_;

  /**
   * Construct a variable implementation from a value. The
   * adjoint is initialized to zero.
   *
   * All constructed variables are added to the stack. For a sparse eigen matrix
   * this includes the nozero values as well the inner and outer indices.
   * Variables should be constructed before variables on which they depend
   * to insure proper partial derivative propagation.  During
   * derivative propagation, the chain() method of each variable
   * will be called in the reverse order of construction.
   *
   * @tparam S A sparse Eigen type that is convertible to `value_type`
   * @param x Value of the constructed variable.
   */
  template <typename S, require_convertible_t<S&, value_type>* = nullptr>
  explicit vari_value(S&& x) : val_(x), adj_(x), chainable_alloc() {
    this->set_zero_adjoint();
    ChainableStack::instance_->var_stack_.push_back(this);
  }
  /**
   * Construct an sparse Eigen variable implementation from a value. The
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
   * @tparam S A sparse Eigen type that is convertible to `value_type`
   * @param x Value of the constructed variable.
   * @param stacked If false will put this this vari on the nochain stack so
   * that its `chain()` method is not called.
   */
  template <typename S, require_convertible_t<S&, value_type>* = nullptr>
  vari_value(S&& x, bool stacked) : val_(x), adj_(x), chainable_alloc() {
    this->set_zero_adjoint();
    if (stacked) {
      ChainableStack::instance_->var_stack_.push_back(this);
    } else {
      ChainableStack::instance_->var_nochain_stack_.push_back(this);
    }
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

  void chain() {}
  /**
   * Initialize the adjoint for this (dependent) variable to 1.
   * This operation is applied to the dependent variable before
   * propagating derivatives, setting the derivative of the
   * result with respect to itself to be 1.
   */
  inline void init_dependent() {
    for (int k = 0; k < adj_.outerSize(); ++k) {
      for (typename PlainObject::InnerIterator it(adj_, k); it; ++it) {
        it.valueRef() = 1.0;
      }
    }
  }

  /**
   * Set the adjoint value of this variable to 0.  This is used to
   * reset adjoints before propagating derivatives again (for
   * example in a Jacobian calculation).
   */
  inline void set_zero_adjoint() noexcept final {
    for (int k = 0; k < adj_.outerSize(); ++k) {
      for (typename PlainObject::InnerIterator it(adj_, k); it; ++it) {
        it.valueRef() = 0.0;
      }
    }
  }

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
    return os << "val: \n" << v->val_ << " \nadj: \n" << v->adj_;
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

 private:
  template <typename>
  friend class var_value;
  template <typename, typename>
  friend class vari_value;
};

}  // namespace math

/**
 * Template specialization defining the scalar type of
 * values stored in vari_value.
 *
 * @tparam T type to check.
 * @ingroup type_trait
 */
template <typename T>
struct scalar_type<math::vari_value<T>> {
  using type = math::vari_value<scalar_type_t<T>>;
};

}  // namespace stan
#endif
