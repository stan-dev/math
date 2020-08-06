#ifndef STAN_MATH_REV_CORE_VAR_HPP
#define STAN_MATH_REV_CORE_VAR_HPP

#include <stan/math/rev/core/vari.hpp>
#include <stan/math/rev/core/grad.hpp>
#include <stan/math/rev/core/chainable_alloc.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/rev/meta/is_vari.hpp>
#include <ostream>
#include <vector>
#ifdef STAN_OPENCL
#include <stan/math/opencl/rev/vari.hpp>
#endif

namespace stan {
namespace math {

// forward declare
template <typename Vari>
static void grad(Vari* vi);

/**
 * Independent (input) and dependent (output) variables for gradients.
 *
 * This class acts as a smart pointer, with resources managed by
 * an arena-based memory manager scoped to a single gradient
 * calculation.
 *
 * A var is constructed with a type `T` and used like any
 * other scalar. Arithmetical functions like negation, addition,
 * and subtraction, as well as a range of mathematical functions
 * like exponentiation and powers are overridden to operate on
 * var values objects.
 * @tparam T An Floating point type.
 */
template <typename T>
class var_value {
  static_assert(
      is_plain_type<T>::value,
      "The template for this var is an"
      " expression but a var_value's inner type must be assignable such as"
      " a double, Eigen::Matrix, or Eigen::Array");
  static_assert(
      std::is_floating_point<value_type_t<T>>::value,
      "The template for must be a floating point or a container holding"
      " floating point types");

 public:
  using value_type = std::decay_t<T>;        // type in vari_value.
  using vari_type = vari_value<value_type>;  // Type of underlying vari impl.

  /**
   * Pointer to the implementation of this variable.
   *
   * This value should not be modified, but may be accessed in
   * <code>var</code> operators to construct `vari_value<T>`
   * instances.
   */
  vari_type* vi_;

  /**
   * Return `true` if this variable has been
   * declared, but not been defined.  Any attempt to use an
   * undefined variable's value or adjoint will result in a
   * segmentation fault.
   *
   * @return <code>true</code> if this variable does not yet have
   * a defined variable.
   */
  inline bool is_uninitialized() { return (vi_ == nullptr); }

  /**
   * Construct a variable for later assignment.
   *
   * This is implemented as a no-op, leaving the underlying implementation
   * dangling.  Before an assignment, the behavior is thus undefined just
   * as for a basic double.
   */
  var_value() : vi_(nullptr) {}

  /**
   * Construct a variable from the specified floating point argument
   * by constructing a new `vari_value<value_type>`. This constructor is only
   * valid when `S` is convertible to this `vari_value`'s `value_type`.
   * @tparam S A type that is convertible to `value_type`.
   * @param x Value of the variable.
   */
  template <typename S, require_convertible_t<S&, value_type>* = nullptr>
  var_value(const S& x) : vi_(new vari_type(x, false)) {}  // NOLINT

  /**
   * Construct a variable from a pointer to a variable implementation.
   * @param vi A vari_value pointer.
   */
  var_value(vari_type* vi) : vi_(vi) {}  // NOLINT

  /**
   * Return a constant reference to the value of this variable.
   *
   * @return The value of this variable.
   */
  inline const auto& val() const { return vi_->val_; }

  /**
   * Return a const reference of the derivative of the root expression with
   * respect to this expression.  This method only works
   * after one of the `grad()` methods has been
   * called.
   *
   * @return Adjoint for this variable.
   */
  inline const auto& adj() const { return vi_->adj_; }

  /**
   * Return a reference to the derivative of the root expression with
   * respect to this expression.  This method only works
   * after one of the `grad()` methods has been
   * called.
   *
   * @return Adjoint for this variable.
   */
  inline auto& adj() { return vi_->adj_; }

  /**
   * Compute the gradient of this (dependent) variable with respect to
   * the specified vector of (independent) variables, assigning the
   * specified vector to the gradient.
   *
   * The grad() function does <i>not</i> recover memory.  In Stan
   * 2.4 and earlier, this function did recover memory.
   *
   * @tparam CheckContainer Not set by user. The default value of value_type
   *  is used to require that grad is only available for scalar `var_value`
   *  types.
   * @param x Vector of independent variables.
   * @param g Gradient vector of partial derivatives of this
   * variable with respect to x.
   */
  template <typename CheckContainer = value_type,
            require_not_container_t<CheckContainer>* = nullptr>
  inline void grad(std::vector<var_value<T>>& x, std::vector<value_type>& g) {
    stan::math::grad(vi_);
    g.resize(x.size());
    for (size_t i = 0; i < x.size(); ++i) {
      g[i] = x[i].vi_->adj_;
    }
  }

  /**
   * Compute the gradient of this (dependent) variable with respect
   * to all (independent) variables.
   *
   * @tparam CheckContainer Not set by user. The default value of value_type
   *  is used to require that grad is only available for scalar `var_value`
   *  types.
   * The grad() function does <i>not</i> recover memory.
   */
  template <typename CheckContainer = value_type,
            require_not_container_t<CheckContainer>* = nullptr>
  void grad() {
    stan::math::grad(vi_);
  }

  // POINTER OVERRIDES

  /**
   * Return a reference to underlying implementation of this variable.
   *
   * If <code>x</code> is of type <code>var</code>, then applying
   * this operator, <code>*x</code>, has the same behavior as
   * <code>*(x.vi_)</code>.
   *
   * <i>Warning</i>:  The returned reference does not track changes to
   * this variable.
   *
   * @return variable
   */
  inline vari_type& operator*() { return *vi_; }

  /**
   * Return a pointer to the underlying implementation of this variable.
   *
   * If <code>x</code> is of type <code>var</code>, then applying
   * this operator, <code>x-&gt;</code>, behaves the same way as
   * <code>x.vi_-&gt;</code>.
   *
   * <i>Warning</i>: The returned result does not track changes to
   * this variable.
   */
  inline vari_type* operator->() { return vi_; }

  // COMPOUND ASSIGNMENT OPERATORS

  /**
   * The compound add/assignment operator for variables (C++).
   *
   * If this variable is a and the argument is the variable b,
   * then (a += b) behaves exactly the same way as (a = a + b),
   * creating an intermediate variable representing (a + b).
   *
   * @param b The variable to add to this variable.
   * @return The result of adding the specified variable to this variable.
   */
  inline var_value<T>& operator+=(const var_value<T>& b);

  /**
   * The compound add/assignment operator for scalars (C++).
   *
   * If this variable is a and the argument is the scalar b, then
   * (a += b) behaves exactly the same way as (a = a + b).  Note
   * that the result is an assignable lvalue.
   *
   * @param b The scalar to add to this variable.
   * @return The result of adding the specified variable to this variable.
   */
  inline var_value<T>& operator+=(T b);

  /**
   * The compound subtract/assignment operator for variables (C++).
   *
   * If this variable is a and the argument is the variable b,
   * then (a -= b) behaves exactly the same way as (a = a - b).
   * Note that the result is an assignable lvalue.
   *
   * @param b The variable to subtract from this variable.
   * @return The result of subtracting the specified variable from
   * this variable.
   */
  inline var_value<T>& operator-=(const var_value<T>& b);

  /**
   * The compound subtract/assignment operator for scalars (C++).
   *
   * If this variable is a and the argument is the scalar b, then
   * (a -= b) behaves exactly the same way as (a = a - b).  Note
   * that the result is an assignable lvalue.
   *
   * @param b The scalar to subtract from this variable.
   * @return The result of subtracting the specified variable from this
   * variable.
   */
  inline var_value<T>& operator-=(T b);

  /**
   * The compound multiply/assignment operator for variables (C++).
   *
   * If this variable is a and the argument is the variable b,
   * then (a *= b) behaves exactly the same way as (a = a * b).
   * Note that the result is an assignable lvalue.
   *
   * @param b The variable to multiply this variable by.
   * @return The result of multiplying this variable by the
   * specified variable.
   */
  inline var_value<T>& operator*=(const var_value<T>& b);

  /**
   * The compound multiply/assignment operator for scalars (C++).
   *
   * If this variable is a and the argument is the scalar b, then
   * (a *= b) behaves exactly the same way as (a = a * b).  Note
   * that the result is an assignable lvalue.
   *
   * @param b The scalar to multiply this variable by.
   * @return The result of multiplying this variable by the specified
   * variable.
   */
  inline var_value<T>& operator*=(T b);

  /**
   * The compound divide/assignment operator for variables (C++).  If this
   * variable is a and the argument is the variable b, then (a /= b)
   * behaves exactly the same way as (a = a / b).  Note that the
   * result is an assignable lvalue.
   *
   * @param b The variable to divide this variable by.
   * @return The result of dividing this variable by the
   * specified variable.
   */
  inline var_value<T>& operator/=(const var_value<T>& b);

  /**
   * The compound divide/assignment operator for scalars (C++).
   *
   * If this variable is a and the argument is the scalar b, then
   * (a /= b) behaves exactly the same way as (a = a / b).  Note
   * that the result is an assignable lvalue.
   *
   * @param b The scalar to divide this variable by.
   * @return The result of dividing this variable by the specified
   * variable.
   */
  inline var_value<T>& operator/=(T b);

  /**
   * Write the value of this autodiff variable and its adjoint to
   * the specified output stream.
   *
   * @param os Output stream to which to write.
   * @param v Variable to write.
   * @return Reference to the specified output stream.
   */
  friend std::ostream& operator<<(std::ostream& os, const var_value<T>& v) {
    if (v.vi_ == nullptr) {
      return os << "uninitialized";
    }
    return os << v.val();
  }
};

// For backwards compatability the default value is double
using var = var_value<double>;

}  // namespace math
}  // namespace stan
#endif
