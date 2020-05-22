#ifndef STAN_MATH_REV_CORE_VAR_HPP
#define STAN_MATH_REV_CORE_VAR_HPP

#include <stan/math/rev/core/vari.hpp>
#include <stan/math/rev/core/grad.hpp>
#include <stan/math/rev/core/chainable_alloc.hpp>
#include <stan/math/rev/core/propogate_static_matrix.hpp>
#include <stan/math/prim/meta.hpp>
#include <boost/math/tools/config.hpp>
#include <ostream>
#include <vector>
#include <complex>
#include <string>
#include <exception>

namespace stan {
namespace math {

// forward declare
static void grad(vari_base* vi);

/**
 * Independent (input) and dependent (output) variables for gradients.
 *
 * This class acts as a smart pointer, with resources managed by
 * an arena-based memory manager scoped to a single gradient
 * calculation.
 *
 * A var is constructed with a double and used like any
 * other scalar.  Arithmetical functions like negation, addition,
 * and subtraction, as well as a range of mathematical functions
 * like exponentiation and powers are overridden to operate on
 * var values objects.
 */
template <typename T>
class var_value {
  template <typename Val>
  using floating_point_promoter
      = std::conditional_t<std::is_integral<std::decay_t<Val>>::value, double,
                           std::decay_t<Val>>;

 public:
  // FIXME: doc what this is for
  using Scalar = floating_point_promoter<T>;
  using vari_type = vari_value<T>;
  /**
   * Pointer to the implementation of this variable.
   *
   * This value should not be modified, but may be accessed in
   * <code>var</code> operators to construct `vari_value<T>`
   * instances.
   */
  vari_value<T>* vi_;

  /**
   * Return `true` if this variable has been
   * declared, but not been defined.  Any attempt to use an
   * undefined variable's value or adjoint will result in a
   * segmentation fault.
   *
   * @return <code>true</code> if this variable does not yet have
   * a defined variable.
   */
  bool is_uninitialized() {
    return (vi_ == static_cast<vari_value<T>*>(nullptr));
  }

  /**
   * Construct a variable for later assignment.
   *
   * This is implemented as a no-op, leaving the underlying implementation
   * dangling.  Before an assignment, the behavior is thus undefined just
   * as for a basic double.
   */
  var_value() : vi_(static_cast<vari_value<T>*>(nullptr)) {}

  /**
   * Construct a variable from a pointer to a variable implementation.
   *
   * @param vi Variable implementation.
   */
  var_value(vari_value<T>* vi) : vi_(vi) {}  // NOLINT

  /**
   * Construct a variable from the specified integral type argument
   * by constructing a new `vari_value<T>`. For integral types the
   * `vari_value<T>` will hold doubles. This constructor is only valid when `T`
   * is arithmetic.
   *
   * @param x Value of the variable.
   */
  template <typename IntegralT, typename T1 = T,
            require_not_same_t<T1, IntegralT>* = nullptr,
            require_integral_t<IntegralT>* = nullptr,
            require_arithmetic_t<T1>* = nullptr>
  var_value(IntegralT x) : vi_(new vari_value<T>(x, false)) {}  // NOLINT

  var_value(T x) : vi_(new vari_value<T>(x, false)) {}  // NOLINT

  template <typename EigenT, typename T1 = T,
            require_not_same_t<T1, EigenT>* = nullptr,
            require_all_eigen_t<EigenT, T1>* = nullptr,
            require_eigen_vt<std::is_arithmetic, EigenT>* = nullptr>
  var_value(EigenT x) : vi_(new vari_value<T>(x, false)) {}  // NOLINT

  template <typename EigenT, typename T1 = T,
            require_not_same_t<T1, EigenT>* = nullptr,
            require_all_eigen_t<EigenT, T1>* = nullptr,
            require_eigen_vt<is_var, EigenT>* = nullptr>
  var_value(EigenT x);  // NOLINT

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
  inline void grad(std::vector<var_value<T>>& x, std::vector<Scalar>& g) {
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
  inline vari_value<T>& operator*() { return *vi_; }

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
  inline vari_value<T>* operator->() { return vi_; }

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
  template <typename Arith, require_arithmetic_t<Arith>...>
  inline var_value<T>& operator+=(const Arith& b);

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
  template <typename Arith, require_arithmetic_t<Arith>...>
  inline var_value<T>& operator-=(const Arith& b);

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
  template <typename Arith, require_vt_arithmetic<Arith>...>
  inline var_value<T>& operator*=(const Arith& b);

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
  template <typename Arith, require_arithmetic_t<Arith>...>
  inline var_value<T>& operator/=(const Arith& b);

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

  template <int R, int C>
  operator Eigen::Matrix<var_value<double>, R, C>();
};

template <typename T>
template <typename EigenT, typename T1, require_not_same_t<T1, EigenT>*,
          require_all_eigen_t<EigenT, T1>*, require_eigen_vt<is_var, EigenT>*>
var_value<T>::var_value(EigenT x)
    : vi_(new vari_value<T>(x.val(), false)) {  // NOLINT
  ChainableStack::instance_->var_stack_.push_back(
      new static_to_dynamic_vari<T>(x.data()[0].vi_, this->vi_, x.size()));
}

template <typename T>
template <int R, int C>
var_value<T>::operator Eigen::Matrix<var_value<double>, R, C>() {
  Eigen::Matrix<var_value<double>, R, C> x(this->val());
  ChainableStack::instance_->var_stack_.push_back(
      new dynamic_to_static_vari<T>(this->vi_, x.data()[0].vi_, x.size()));
  return x;
}

using var = var_value<double>;

}  // namespace math
}  // namespace stan
#endif
