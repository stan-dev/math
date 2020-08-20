#ifndef STAN_MATH_REV_CORE_VARI_HPP
#define STAN_MATH_REV_CORE_VARI_HPP

#include <stan/math/rev/core/chainable_alloc.hpp>
#include <stan/math/rev/core/chainablestack.hpp>
#include <stan/math/prim/meta.hpp>
#include <ostream>
#include <type_traits>

namespace stan {
namespace math {

// forward decleration of vari_value
template <typename T, typename StrideType = Eigen::Stride<0, 0>, typename = void>
class vari_value;

// forward declaration of var
template <typename T,
          typename VariType = vari_value<std::remove_reference_t<T>>>
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
template <typename T, typename StrideType>
class vari_value<T, StrideType, require_floating_point_t<T>> : public vari_base {
 public:
  using value_type = std::decay_t<T>;
  /**
   * The value of this variable.
   */
  const value_type val_;
  /**
   * The adjoint of this variable, which is the partial derivative
   * of this variable with respect to the root variable.
   */
  value_type adj_;

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
  template <typename S, require_convertible_t<S&, T>* = nullptr>
  vari_value(S x) noexcept : val_(x), adj_(0.0) {  // NOLINT
    ChainableStack::instance_->var_stack_.push_back(this);
  }

  /**
   * Construct a variable implementation from a value.  The
   *  adjoint is initialized to zero and if `stacked` is `false` this vari
   *  will be not be put on the var_stack. Instead it will only be put on
   *  a stack to keep track of whether the adjoint needs to be set to zero.
   *  Variables should be constructed before variables on which they depend
   *  to insure proper partial derivative propagation.  During
   *  derivative propagation, the chain() method of each variable
   *  will be called in the reverse order of construction.
   *
   * @tparam S n floating point type.
   * @param x Value of the constructed variable.
   * @param stacked If false will put this this vari on the nochain stack so
   * that its `chain()` method is not called.
   */
  template <typename S, require_convertible_t<S&, T>* = nullptr>
  vari_value(S x, bool stacked) noexcept : val_(x), adj_(0.0) {
    if (stacked) {
      ChainableStack::instance_->var_stack_.push_back(this);
    } else {
      ChainableStack::instance_->var_nochain_stack_.push_back(this);
    }
  }

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
  friend std::ostream& operator<<(std::ostream& os, const vari_value<T, StrideType>* v) {
    return os << v->val_ << ":" << v->adj_;
  }
protected:
  template <typename S, typename K, require_all_floating_point_t<S, K>* = nullptr>
  vari_value(S val, K adj) : val_(val), adj_(adj) {}

 private:
  template <typename, typename>
  friend class var_value;
  template <typename, typename, typename>
  friend class vari_value;

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
template <typename T, typename StrideType>
class vari_value<T, StrideType, require_eigen_dense_base_t<T>> : public vari_base {
 public:

  static_assert(
      is_plain_type<T>::value,
      "The template for this var is an"
      " expression but a var_value's inner type must be assignable such as"
      " a double, Eigen::Matrix, or Eigen::Array");

  /**
   * `PlainObject` represents a user constructible type such as Matrix or Array
   */
  using PlainObject = plain_type_t<T>;
  using value_type = PlainObject;  // The underlying type for this class
  using eigen_scalar = value_type_t<PlainObject>;  // A floating point type
  using eigen_map
    = Eigen::Map<T, Eigen::Aligned8, StrideType>;  // Maps for adj_ and val_
  using Stride = StrideType;
  using vari_type = vari_value<T, StrideType>;
  /**
   * Number of rows known at compile time
   */
  static constexpr int RowsAtCompileTime = PlainObject::RowsAtCompileTime;
  /**
   * Number of columns known at compile time
   */
  static constexpr int ColsAtCompileTime = PlainObject::ColsAtCompileTime;
  /**
   * Maps for adj_ and val_
   */
  eigen_scalar* val_mem_;  // Pointer to memory allocated on the stack for val_
  eigen_scalar* adj_mem_;  // Pointer to memory allocated on the stack for adj_

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
  template <typename S, require_convertible_t<S&, T>* = nullptr>
  explicit vari_value(const S& x)
      : val_mem_(ChainableStack::instance_->memalloc_.alloc_array<eigen_scalar>(
            x.size())),
        adj_mem_(ChainableStack::instance_->memalloc_.alloc_array<eigen_scalar>(
            x.size())),
        val_(eigen_map(val_mem_, x.rows(), x.cols()) = x),
        adj_(eigen_map(adj_mem_, x.rows(), x.cols()).setZero()) {
    ChainableStack::instance_->var_stack_.push_back(this);
  }

  /**
   * Construct a dense Eigen variable implementation from a value. The
   *  adjoint is initialized to zero and if `stacked` is `false` this vari
   *  will be not be put on the var_stack. Instead it will only be put on
   *  a stack to keep track of whether the adjoint needs to be set to zero.
   *  Variables should be constructed before variables on which they depend
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
  vari_value(const S& x, bool stacked)
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
   * Implimentation for assignable objects that sets the adjoint values to zero.
   */
  template <typename Checker = T, require_t<bool_constant<!std::is_const<Checker>::value>>* = nullptr>
  inline void set_zero_adjoint_impl() {
    adj_.setZero();
  }
  /**
   * Implimentation for views which performs a no-op
   */
  template <typename Checker = T, require_t<std::is_const<Checker>>* = nullptr>
  inline void set_zero_adjoint_impl() {}

  /**
   * Set the adjoint value of this variable to 0.  This is used to
   * reset adjoints before propagating derivatives again (for
   * example in a Jacobian calculation).
   */
  inline void set_zero_adjoint() final { set_zero_adjoint_impl(); }

    /**
     * A block view of the underlying Eigen matrices.
     * @param i Starting row of block.
     * @param j Starting columns of block.
     * @param p Number of rows to return.
     * @param q Number of columns to return.
     */
  inline const auto block(Eigen::Index i, Eigen::Index j, Eigen::Index p,
                     Eigen::Index q) const {
       const auto& val_block = val_.block(i, j, p, q);
       const auto& adj_block = adj_.block(i, j, p, q);
      return vari_value<const PlainObject, Eigen::OuterStride<>>(
        val_block, adj_block, val_block.outerStride(), adj_block.outerStride());
    }

    /**
     * View of the head of Eigen vector types.
     * @param n Number of elements to return from top of vector.
     */
    inline const auto head(Eigen::Index n) const {
      return vari_value<const PlainObject>(val_.head(n), adj_.head(n));
    }

    /**
     * View of the tail of the Eigen vector types.
     * @param n Number of elements to return from bottom of vector.
     */
    inline const auto tail(Eigen::Index n) const {
      return vari_value<const PlainObject>(val_.tail(n), adj_.tail(n));
    }

    /**
     * View block of N elements starting at position `i`
     * @param i Starting position of block.
     * @param n Number of elements in block
     */
    const auto segment(Eigen::Index i, Eigen::Index n) const {
      return vari_value<const PlainObject>(val_.segment(i, n), adj_.segment(i, n));
    }

    /**
     * View row of eigen matrices.
     * @param i Row index to slice.
     */
    inline const auto row(Eigen::Index i) const {
      const auto& val_row = val_.row(i);
      const auto& adj_row = adj_.row(i);
      return vari_value<const PlainObject, Eigen::InnerStride<>>(
          val_row, adj_row, val_row.innerStride(), adj_row.innerStride());
    }

    /**
     * View column of eigen matrices
     * @param i Column index to slice
     */
    inline const auto col(Eigen::Index i) const {
      return vari_value<const PlainObject, Eigen::OuterStride<>>(val_.col(i), adj_.col(i),
    val_.col(i).outerStride(), adj_.col(i).outerStride());
    }

    /**
     * Get coefficient of eigen matrices
     * @param i Row index
     * @param j Column index
     */
    inline const auto coeff(Eigen::Index i, Eigen::Index j) const {
      return vari_value<double>(val_(i, j), adj_(i, j));
    }

    /**
     * Get coefficient of eigen matrices
     * @param i Column index to slice
     */
    inline const auto operator()(Eigen::Index i) const {
      return vari_value<double>(val_(i), adj_(i));
    }

    /**
     * Get coefficient of eigen matrices
     * @param i Row index
     * @param j Column index
     */
    inline const auto operator()(Eigen::Index i, Eigen::Index j) const {
      return this->coeff(i, j);
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
  friend std::ostream& operator<<(std::ostream& os,
     const vari_value<T, StrideType>* v) {
    return os << "val: \n" << v->val_ << " \nadj: \n" << v->adj_;
  }
 protected:
   template <typename S, typename K, typename StrideS = Eigen::Stride<0, 0>,
    typename StrideK = Eigen::Stride<0, 0>,
             require_convertible_t<S&, value_type>* = nullptr,
             require_convertible_t<K&, value_type>* = nullptr>
   explicit vari_value(S&& val, K&& adj,
      StrideS val_stride = Eigen::Stride<0, 0>(),
      StrideK adj_stride = Eigen::Stride<0, 0>())
       : val_(val.data(), val.rows(), val.cols(), val_stride),
         adj_(adj.data(), adj.rows(), adj.cols(), adj_stride) {}

 private:
  template <typename, typename>
  friend class var_value;
  template <typename, typename, typename>
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
template <typename T, typename StrideType>
class vari_value<T, StrideType, require_eigen_sparse_base_t<T>> : public vari_base,
                                                      chainable_alloc {
 public:
  using PlainObject = plain_type_t<T>;  // Base type of Eigen class
  using value_type = PlainObject;       // vari's adj_ and val_ member type
  /**
   * Rows at compile time
   */
  static constexpr int RowsAtCompileTime = T::RowsAtCompileTime;
  /**
   * Columns at compile time
   */
  static constexpr int ColsAtCompileTime = T::ColsAtCompileTime;

  /**
   * The adjoint of this variable, which is the partial derivative
   * of this variable with respect to the root variable.
   */
  PlainObject adj_;

  /**
   * The value of this variable.
   */
  const PlainObject val_;

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
  template <typename S, require_convertible_t<S&, T>* = nullptr>
  explicit vari_value(S&& x)
      : adj_(x), val_(std::forward<S>(x)), chainable_alloc() {
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
  template <typename S, require_convertible_t<S&, T>* = nullptr>
  vari_value(S&& x, bool stacked)
      : adj_(x), val_(std::forward<S>(x)), chainable_alloc() {
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
  friend std::ostream& operator<<(std::ostream& os, const vari_value<T, StrideType>* v) {
    return os << "val: \n" << v->val_ << " \nadj: \n" << v->adj_;
  }

 private:
  template <typename, typename>
  friend class var_value;
  template <typename, typename, typename>
  friend class vari_value;
};

// For backwards compatability the default is double
using vari = vari_value<double>;

}  // namespace math
}  // namespace stan
#endif
