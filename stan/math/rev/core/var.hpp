#ifndef STAN_MATH_REV_CORE_VAR_HPP
#define STAN_MATH_REV_CORE_VAR_HPP

#include <stan/math/rev/core/vari.hpp>
#include <stan/math/rev/core/grad.hpp>
#include <stan/math/rev/core/chainable_alloc.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/rev/meta/is_vari.hpp>
#include <ostream>
#include <vector>

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
 * other scalar.  Arithmetical functions like negation, addition,
 * and subtraction, as well as a range of mathematical functions
 * like exponentiation and powers are overridden to operate on
 * var values objects.
 * @tparam T An Arithmetic type.
 */
template <typename T>
class var_value {
 public:
  // The internal value_type is always a floating point type
  using value_type = internal::floating_point_promoter<T>;
  using vari_type = vari_value<value_type>;

  template <typename K>
  using require_vari_convertible_t
      = require_t<std::is_convertible<vari_value<K>*, vari_type*>>;
  template <typename K>
  using require_not_vari_convertible_t = require_t<
      bool_constant<!std::is_convertible<vari_value<K>*, vari_type*>::value>>;

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
  bool is_uninitialized() { return (vi_ == static_cast<vari_type*>(nullptr)); }

  /**
   * Construct a variable for later assignment.
   *
   * This is implemented as a no-op, leaving the underlying implementation
   * dangling.  Before an assignment, the behavior is thus undefined just
   * as for a basic double.
   */
  var_value() : vi_(static_cast<vari_type*>(nullptr)) {}

  /**
   * Construct a variable from the specified arithmetic type argument
   * by constructing a new `vari_value<value_type>`. For integral types the
   * `vari_value<value_type>` will hold doubles. This constructor is only valid
   * when `S` is convertible to this `vari_value`'s `value_type`.
   * @tparam S A type that is convertible to `value_type`.
   * @param x Value of the variable.
   */
  template <typename S, require_convertible_t<S&, value_type>* = nullptr>
  var_value(S x) : vi_(new vari_type(x, false)) {}  // NOLINT

  /**
   * Construct a variable from a pointer to a variable implementation.
   * @tparam S The type in the vari_value pointer that has the same `value_type`
   *  as this `var_value`. For integral types a `vari_value<Integral>` is the
   * same as `vari_value<double>` so can be `reinterpret_cast` without a copy.
   * @param vi A vari_value pointer.
   */
  template <typename S, typename K = internal::floating_point_promoter<S>,
            require_vari_convertible_t<K>* = nullptr>
  var_value(vari_value<S>* vi) // NOLINT
      : vi_(reinterpret_cast<vari_type*>(vi)) {}

  /**
   * Construct a variable from a pointer to a variable implementation when.
   *  the `value_type`s differ. This forces a copy of the vari
   * @tparam S the type in the vari_value* that does not have the same type
   *  as the `value_type` of this `var_value`.
   * @param vi A `vari_value`.
   */
  template <typename S, typename K = internal::floating_point_promoter<S>,
            require_not_vari_convertible_t<K>* = nullptr>
  var_value(vari_value<S>* vi) : vi_(new vari_type(*vi)) {}  // NOLINT

  /**
   * Construct from a var_value whose underlying value_type differs in type
   * to this class's value_type. This is for cases such as
   * `var_value<double> a(5.0); var_value<long double> b(a);`
   * @tparam S A type that is not the same as the underlying value type.
   *  same as the underlying value_type.
   * param x a `var_value` whose underlying vari_type can be dynamically cast
   * to `this::vari_value<value_type>``.
   */
  template <typename S, typename K = internal::floating_point_promoter<S>,
            require_not_vari_convertible_t<K>* = nullptr>
  var_value(var_value<S>& x) : vi_(new vari_type(*x.vi_)) {}  // NOLINT

  template <typename S, typename K = internal::floating_point_promoter<S>,
            require_not_vari_convertible_t<K>* = nullptr>
  var_value(var_value<S>&& x) : vi_(new vari_type(*x.vi_)) {}  // NOLINT

  /**
   * Constructor from `var_value` whose value_type is the same as this class's
   * `value_type`. This is used in cases such as
   * `var_value<double> a(4.0); var_value<int> b(a)` since the `value_type` for
   *  a `var_value` with an integral type is a double.
   */
  template <typename S, typename K = internal::floating_point_promoter<S>,
            require_vari_convertible_t<K>* = nullptr>
  var_value(var_value<S>& x) : vi_(x.vi_) {}  // NOLINT
  template <typename S, typename K = internal::floating_point_promoter<S>,
            require_vari_convertible_t<K>* = nullptr>
  var_value(var_value<S>&& x) : vi_(x.vi_) {}  // NOLINT

  /**
   * Return the value of this variable.
   *
   * @return The value of this variable.
   */
  inline auto val() const { return vi_->val_; }

  /**
   * Return the derivative of the root expression with
   * respect to this expression.  This method only works
   * after one of the `grad()` methods has been
   * called.
   *
   * @return Adjoint for this variable.
   */
  inline auto adj() const { return vi_->adj_; }

  /**
   * Compute the gradient of this (dependent) variable with respect to
   * the specified vector of (independent) variables, assigning the
   * specified vector to the gradient.
   *
   * The grad() function does <i>not</i> recover memory.  In Stan
   * 2.4 and earlier, this function did recover memory.
   *
   * @param x Vector of independent variables.
   * @param g Gradient vector of partial derivatives of this
   * variable with respect to x.
   */
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
   * The grad() function does <i>not</i> recover memory.
   */
  void grad() { stan::math::grad(vi_); }

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
   * @tparam S a type that is convertible to `value_type`
   * @param b The variable to add to this variable.
   * @return The result of adding the specified variable to this variable.
   */
  template <typename S,
            require_convertible_t<
                S&, internal::floating_point_promoter<T>>* = nullptr>
  inline var_value<T>& operator+=(const var_value<S>& b);

  /**
   * The compound add/assignment operator for scalars (C++).
   *
   * If this variable is a and the argument is the scalar b, then
   * (a += b) behaves exactly the same way as (a = a + b).  Note
   * that the result is an assignable lvalue.
   *
   * @tparam Arith An arithmetic type
   * @param b The scalar to add to this variable.
   * @return The result of adding the specified variable to this variable.
   */
  template <typename Arith, require_arithmetic_t<Arith>* = nullptr>
  inline var_value<T>& operator+=(Arith b);

  /**
   * The compound subtract/assignment operator for variables (C++).
   *
   * If this variable is a and the argument is the variable b,
   * then (a -= b) behaves exactly the same way as (a = a - b).
   * Note that the result is an assignable lvalue.
   *
   * @tparam S a type that is convertible to `value_type`
   * @param b The variable to subtract from this variable.
   * @return The result of subtracting the specified variable from
   * this variable.
   */
  template <typename S,
            require_convertible_t<
                S&, internal::floating_point_promoter<T>>* = nullptr>
  inline var_value<T>& operator-=(const var_value<S>& b);

  /**
   * The compound subtract/assignment operator for scalars (C++).
   *
   * If this variable is a and the argument is the scalar b, then
   * (a -= b) behaves exactly the same way as (a = a - b).  Note
   * that the result is an assignable lvalue.
   *
   * @tparam Arith An arithmetic type
   * @param b The scalar to subtract from this variable.
   * @return The result of subtracting the specified variable from this
   * variable.
   */
  template <typename Arith, require_arithmetic_t<Arith>* = nullptr>
  inline var_value<T>& operator-=(Arith b);

  /**
   * The compound multiply/assignment operator for variables (C++).
   *
   * If this variable is a and the argument is the variable b,
   * then (a *= b) behaves exactly the same way as (a = a * b).
   * Note that the result is an assignable lvalue.
   *
   * @tparam S a type that is convertible to `value_type`
   * @param b The variable to multiply this variable by.
   * @return The result of multiplying this variable by the
   * specified variable.
   */
  template <typename S,
            require_convertible_t<
                S&, internal::floating_point_promoter<T>>* = nullptr>
  inline var_value<T>& operator*=(const var_value<S>& b);

  /**
   * The compound multiply/assignment operator for scalars (C++).
   *
   * If this variable is a and the argument is the scalar b, then
   * (a *= b) behaves exactly the same way as (a = a * b).  Note
   * that the result is an assignable lvalue.
   *
   * @tparam Arith An arithmetic type
   * @param b The scalar to multiply this variable by.
   * @return The result of multiplying this variable by the specified
   * variable.
   */
  template <typename Arith, require_arithmetic_t<Arith>* = nullptr>
  inline var_value<T>& operator*=(Arith b);

  /**
   * The compound divide/assignment operator for variables (C++).  If this
   * variable is a and the argument is the variable b, then (a /= b)
   * behaves exactly the same way as (a = a / b).  Note that the
   * result is an assignable lvalue.
   *
   * @tparam S a type that is convertible to `value_type`
   * @param b The variable to divide this variable by.
   * @return The result of dividing this variable by the
   * specified variable.
   */
  template <typename S,
            require_convertible_t<
                S&, internal::floating_point_promoter<T>>* = nullptr>
  inline var_value<T>& operator/=(const var_value<S>& b);

  /**
   * The compound divide/assignment operator for scalars (C++).
   *
   * If this variable is a and the argument is the scalar b, then
   * (a /= b) behaves exactly the same way as (a = a / b).  Note
   * that the result is an assignable lvalue.
   *
   * @tparam Arith An arithmetic type
   * @param b The scalar to divide this variable by.
   * @return The result of dividing this variable by the specified
   * variable.
   */
  template <typename Arith, require_arithmetic_t<Arith>* = nullptr>
  inline var_value<T>& operator/=(Arith b);

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
