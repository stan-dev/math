#ifndef STAN_MATH_REV_CORE_VAR_HPP
#define STAN_MATH_REV_CORE_VAR_HPP

#ifdef STAN_OPENCL
#include <stan/math/opencl/rev/vari.hpp>
#include <stan/math/opencl/plain_type.hpp>
#endif
#include <stan/math/rev/core/vari.hpp>
#include <stan/math/rev/core/grad.hpp>
#include <stan/math/rev/core/chainable_alloc.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/rev/meta/is_vari.hpp>
#include <stan/math/rev/meta/arena_type.hpp>
#include <stan/math/rev/core/reverse_pass_callback.hpp>
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
 * other scalar. Arithmetical functions like negation, addition,
 * and subtraction, as well as a range of mathematical functions
 * like exponentiation and powers are overridden to operate on
 * var values objects.
 * @tparam T An Floating point type.
 */
template <typename T>
class var_value<T, require_floating_point_t<T>> {
 public:
  using value_type = std::decay_t<T>;  // type in vari_value.
  using vari_type = vari_value<value_type>;

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
  var_value(S x) : vi_(new vari_type(x, false)) {}  // NOLINT

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
  inline const auto& val() const noexcept { return vi_->val(); }

  /**
   * Return a reference of the derivative of the root expression with
   * respect to this expression.  This method only works
   * after one of the `grad()` methods has been
   * called.
   *
   * @return Adjoint for this variable.
   */
  inline auto& adj() const noexcept { return vi_->adj(); }

  /**
   * Return a reference to the derivative of the root expression with
   * respect to this expression.  This method only works
   * after one of the `grad()` methods has been
   * called.
   *
   * @return Adjoint for this variable.
   */
  inline auto& adj() noexcept { return vi_->adj_; }

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

namespace internal {
template <typename T>
using require_matrix_var_value = require_t<bool_constant<
    (is_eigen<T>::value || is_kernel_expression_and_not_scalar<T>::value)
    && std::is_floating_point<value_type_t<T>>::value>>;
}

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
class var_value<T, internal::require_matrix_var_value<T>> {
 public:
  using value_type = T;  // type in vari_value.
  using vari_type = std::conditional_t<is_plain_type<value_type>::value,
                                       vari_value<value_type>, vari_view<T>>;

  static constexpr int RowsAtCompileTime{vari_type::RowsAtCompileTime};
  static constexpr int ColsAtCompileTime{vari_type::ColsAtCompileTime};

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
  inline bool is_uninitialized() noexcept { return (vi_ == nullptr); }

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
  template <typename S, require_assignable_t<value_type, S>* = nullptr>
  var_value(S&& x) : vi_(new vari_type(std::forward<S>(x), false)) {}  // NOLINT

  /**
   * Copy constructor for var_val.
   * @tparam S type of the value in the `var_value` to assing
   * @param other the value to assign
   * @return this
   */
  template <typename S, require_assignable_t<value_type, S>* = nullptr,
            require_all_plain_type_t<T, S>* = nullptr>
  var_value(const var_value<S>& other) : vi_(other.vi_) {}

  /**
   * Construct a `var_value` with a plain type
   *  from another `var_value` containing an expression.
   * @tparam S type of the value in the `var_value` to assing
   * @param other the value to assign
   * @return this
   */
  template <typename S, typename T_ = T,
            require_assignable_t<value_type, S>* = nullptr,
            require_not_plain_type_t<S>* = nullptr,
            require_plain_type_t<T_>* = nullptr>
  var_value(const var_value<S>& other) : vi_(new vari_type(other.vi_->val_)) {
    reverse_pass_callback(
        [this_vi = this->vi_, other_vi = other.vi_]() mutable {
          other_vi->adj_ += this_vi->adj_;
        });
  }

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
  inline const auto& val() const noexcept { return vi_->val(); }
  inline auto& val_op() noexcept { return vi_->val_op(); }

  /**
   * Return a reference to the derivative of the root expression with
   * respect to this expression.  This method only works
   * after one of the `grad()` methods has been
   * called.
   *
   * @return Adjoint for this variable.
   */
  inline auto& adj() noexcept { return vi_->adj(); }
  inline auto& adj() const noexcept { return vi_->adj(); }
  inline auto& adj_op() noexcept { return vi_->adj(); }

  inline Eigen::Index rows() const noexcept { return vi_->rows(); }
  inline Eigen::Index cols() const noexcept { return vi_->cols(); }
  inline Eigen::Index size() const noexcept { return vi_->size(); }

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
  template <typename S, require_st_var<S>* = nullptr>
  inline var_value<T>& operator-=(const S& b);

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
  template <typename S, require_st_arithmetic<S>* = nullptr>
  inline var_value<T>& operator-=(const S& b);

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
   * A block view of the underlying Eigen matrices.
   * @param start_row Starting row of block.
   * @param start_col Starting columns of block.
   * @param num_rows Number of rows to return.
   * @param num_cols Number of columns to return.
   */
  inline auto block(Eigen::Index start_row, Eigen::Index start_col,
                    Eigen::Index num_rows, Eigen::Index num_cols) const {
    using vari_sub
        = decltype(vi_->block(start_row, start_col, num_rows, num_cols));
    using var_sub = var_value<value_type_t<vari_sub>>;
    return var_sub(
        new vari_sub(vi_->block(start_row, start_col, num_rows, num_cols)));
  }
  inline auto block(Eigen::Index start_row, Eigen::Index start_col,
                    Eigen::Index num_rows, Eigen::Index num_cols) {
    using vari_sub
        = decltype(vi_->block(start_row, start_col, num_rows, num_cols));
    using var_sub = var_value<value_type_t<vari_sub>>;
    return var_sub(
        new vari_sub(vi_->block(start_row, start_col, num_rows, num_cols)));
  }

  /**
   * View transpose of eigen matrix.
   */
  inline auto transpose() const {
    using vari_sub = decltype(vi_->transpose());
    using var_sub = var_value<value_type_t<vari_sub>>;
    return var_sub(new vari_sub(vi_->transpose()));
  }
  inline auto transpose() {
    using vari_sub = decltype(vi_->transpose());
    using var_sub = var_value<value_type_t<vari_sub>>;
    return var_sub(new vari_sub(vi_->transpose()));
  }

  /**
   * View of the head of Eigen vector types.
   * @param n Number of elements to return from top of vector.
   */
  inline auto head(Eigen::Index n) const {
    using vari_sub = decltype(vi_->head(n));
    using var_sub = var_value<value_type_t<vari_sub>>;
    return var_sub(new vari_sub(vi_->head(n)));
  }
  inline auto head(Eigen::Index n) {
    using vari_sub = decltype(vi_->head(n));
    using var_sub = var_value<value_type_t<vari_sub>>;
    return var_sub(new vari_sub(vi_->head(n)));
  }

  /**
   * View of the tail of the Eigen vector types.
   * @param n Number of elements to return from bottom of vector.
   */
  inline auto tail(Eigen::Index n) const {
    using vari_sub = decltype(vi_->tail(n));
    using var_sub = var_value<value_type_t<vari_sub>>;
    return var_sub(new vari_sub(vi_->tail(n)));
  }
  inline auto tail(Eigen::Index n) {
    using vari_sub = decltype(vi_->tail(n));
    using var_sub = var_value<value_type_t<vari_sub>>;
    return var_sub(new vari_sub(vi_->tail(n)));
  }

  /**
   * View block of N elements starting at position `i`
   * @param i Starting position of block.
   * @param n Number of elements in block
   */
  inline auto segment(Eigen::Index i, Eigen::Index n) const {
    using vari_sub = decltype(vi_->segment(i, n));
    using var_sub = var_value<value_type_t<vari_sub>>;
    return var_sub(new vari_sub(vi_->segment(i, n)));
  }
  inline auto segment(Eigen::Index i, Eigen::Index n) {
    using vari_sub = decltype(vi_->segment(i, n));
    using var_sub = var_value<value_type_t<vari_sub>>;
    return var_sub(new vari_sub(vi_->segment(i, n)));
  }

  /**
   * View row of eigen matrices.
   * @param i Row index to slice.
   */
  inline auto row(Eigen::Index i) const {
    using vari_sub = decltype(vi_->row(i));
    using var_sub = var_value<value_type_t<vari_sub>>;
    return var_sub(new vari_sub(vi_->row(i)));
  }
  inline auto row(Eigen::Index i) {
    using vari_sub = decltype(vi_->row(i));
    using var_sub = var_value<value_type_t<vari_sub>>;
    return var_sub(new vari_sub(vi_->row(i)));
  }

  /**
   * View column of eigen matrices
   * @param i Column index to slice
   */
  inline auto col(Eigen::Index i) const {
    using vari_sub = decltype(vi_->col(i));
    using var_sub = var_value<value_type_t<vari_sub>>;
    return var_sub(new vari_sub(vi_->col(i)));
  }
  inline auto col(Eigen::Index i) {
    using vari_sub = decltype(vi_->col(i));
    using var_sub = var_value<value_type_t<vari_sub>>;
    return var_sub(new vari_sub(vi_->col(i)));
  }

  /**
   * View diagonal of eigen matrices
   * @param i Column index to slice
   */
  inline auto diagonal() const {
    using vari_sub = decltype(vi_->diagonal());
    using var_sub = var_value<value_type_t<vari_sub>>;
    return var_sub(new vari_sub(vi_->diagonal()));
  }
  inline auto diagonal() {
    using vari_sub = decltype(vi_->diagonal());
    using var_sub = var_value<value_type_t<vari_sub>>;
    return var_sub(new vari_sub(vi_->diagonal()));
  }

  /**
   * View a `matrix_cl` as a column vector.
   */
  inline auto as_column_vector_or_scalar() const {
    using vari_sub = decltype(vi_->as_column_vector_or_scalar());
    using var_sub = var_value<value_type_t<vari_sub>>;
    return var_sub(new vari_sub(vi_->as_column_vector_or_scalar()));
  }
  inline auto as_column_vector_or_scalar() {
    using vari_sub = decltype(vi_->as_column_vector_or_scalar());
    using var_sub = var_value<value_type_t<vari_sub>>;
    return var_sub(new vari_sub(vi_->as_column_vector_or_scalar()));
  }

  /**
   * View element of eigen matrices. This creates a new
   * vari_value<double> so unlike the other views this subset will not
   * have the same adjoints as the original matrix and must be propogated
   * back.
   * @param i Element to access
   */
  inline auto coeff(Eigen::Index i) const {
    using vari_sub = decltype(vi_->coeff(i));
    vari_sub* vari_coeff = new vari_sub(vi_->coeff(i));
    reverse_pass_callback([this_vi = this->vi_, vari_coeff, i]() {
      this_vi->adj_.coeffRef(i) += vari_coeff->adj_;
    });
    return var_value<value_type_t<vari_sub>>(vari_coeff);
  }
  inline auto coeff(Eigen::Index i) {
    using vari_sub = decltype(vi_->coeff(i));
    vari_sub* vari_coeff = new vari_sub(vi_->coeff(i));
    reverse_pass_callback([this_vi = this->vi_, vari_coeff, i]() {
      this_vi->adj_.coeffRef(i) += vari_coeff->adj_;
    });
    return var_value<value_type_t<vari_sub>>(vari_coeff);
  }

  /**
   * View element of eigen matrices. This creates a new
   * vari_value<double> so unlike the other views this subset will not
   * have the same adjoints as the original matrix and must be propogated
   * back.
   * @param i Row to access
   * @param j Column to access
   */
  inline auto coeff(Eigen::Index i, Eigen::Index j) const {
    using vari_sub = decltype(vi_->coeff(i, j));
    vari_sub* vari_coeff = new vari_sub(vi_->coeff(i, j));
    reverse_pass_callback([this_vi = this->vi_, vari_coeff, i, j]() {
      this_vi->adj_.coeffRef(i, j) += vari_coeff->adj_;
    });
    return var_value<value_type_t<vari_sub>>(vari_coeff);
  }
  inline auto coeff(Eigen::Index i, Eigen::Index j) {
    using vari_sub = decltype(vi_->coeff(i, j));
    vari_sub* vari_coeff = new vari_sub(vi_->coeff(i, j));
    reverse_pass_callback([this_vi = this->vi_, vari_coeff, i, j]() {
      this_vi->adj_.coeffRef(i, j) += vari_coeff->adj_;
    });
    return var_value<value_type_t<vari_sub>>(vari_coeff);
  }

  /**
   * View element of eigen matrices. This creates a new
   * vari_value<double> so unlike the other views this subset will not
   * have the same adjoints as the original matrix and must be propogated
   * back.
   * @param i Element to access
   */
  inline auto operator()(Eigen::Index i) const { return this->coeff(i); }
  inline auto operator()(Eigen::Index i) { return this->coeff(i); }

  /**
   * View element of eigen matrices. This creates a new
   * vari_value<double> so unlike the other views this subset will not
   * have the same adjoints as the original matrix and must be propogated
   * back.
   * @param i Row to access
   * @param j Column to access
   */
  inline auto operator()(Eigen::Index i, Eigen::Index j) const {
    return this->coeff(i, j);
  }
  inline auto operator()(Eigen::Index i, Eigen::Index j) {
    return this->coeff(i, j);
  }

  /**
   * View element of eigen matrices. This creates a new
   * vari_value<double> so unlike the other views this subset will not
   * have the same adjoints as the original matrix and must be propogated
   * back.
   * @param i Element to access
   */
  inline auto coeffRef(Eigen::Index i) const { return this->coeff(i); }
  inline auto coeffRef(Eigen::Index i) { return this->coeff(i); }

  /**
   * View element of eigen matrices. This creates a new
   * vari_value<double> so unlike the other views this subset will not
   * have the same adjoints as the original matrix and must be propogated
   * back.
   * @param i Row to access
   * @param j Column to access
   */
  inline auto coeffRef(Eigen::Index i, Eigen::Index j) const {
    return this->coeff(i, j);
  }
  inline auto coeffRef(Eigen::Index i, Eigen::Index j) {
    return this->coeff(i, j);
  }

  /**
   * Return an expression that operates on the rows of the matrix `vari`
   */
  inline auto rowwise_reverse() const {
    using vari_sub = decltype(vi_->rowwise_reverse());
    using var_sub = var_value<value_type_t<vari_sub>>;
    return var_sub(new vari_sub(vi_->rowwise_reverse()));
  }
  inline auto rowwise_reverse() {
    using vari_sub = decltype(vi_->rowwise_reverse());
    using var_sub = var_value<value_type_t<vari_sub>>;
    return var_sub(new vari_sub(vi_->rowwise_reverse()));
  }

  /**
   * Return an expression that operates on the columns of the matrix `vari`
   */
  inline auto colwise_reverse() const {
    using vari_sub = decltype(vi_->colwise_reverse());
    using var_sub = var_value<value_type_t<vari_sub>>;
    return var_sub(new vari_sub(vi_->colwise_reverse()));
  }
  inline auto colwise_reverse() {
    using vari_sub = decltype(vi_->colwise_reverse());
    using var_sub = var_value<value_type_t<vari_sub>>;
    return var_sub(new vari_sub(vi_->colwise_reverse()));
  }

  /**
   * Return an expression an expression to reverse the order of the coefficients
   * inside of a `vari` matrix
   */
  inline auto reverse() const {
    using vari_sub = decltype(vi_->reverse());
    using var_sub = var_value<value_type_t<vari_sub>>;
    return var_sub(new vari_sub(vi_->reverse()));
  }
  inline auto reverse() {
    using vari_sub = decltype(vi_->reverse());
    using var_sub = var_value<value_type_t<vari_sub>>;
    return var_sub(new vari_sub(vi_->reverse()));
  }

  /**
   * Return a block consisting of the top rows
   * @param n Number of rows
   */
  inline auto topRows(Eigen::Index n) const {
    using vari_sub = decltype(vi_->topRows(n));
    using var_sub = var_value<value_type_t<vari_sub>>;
    return var_sub(new vari_sub(vi_->topRows(n)));
  }
  inline auto topRows(Eigen::Index n) {
    using vari_sub = decltype(vi_->topRows(n));
    using var_sub = var_value<value_type_t<vari_sub>>;
    return var_sub(new vari_sub(vi_->topRows(n)));
  }

  /**
   * Return a block consisting of the bottom rows
   * @param n Number of rows
   */
  inline auto bottomRows(Eigen::Index n) const {
    using vari_sub = decltype(vi_->bottomRows(n));
    using var_sub = var_value<value_type_t<vari_sub>>;
    return var_sub(new vari_sub(vi_->bottomRows(n)));
  }
  inline auto bottomRows(Eigen::Index n) {
    using vari_sub = decltype(vi_->bottomRows(n));
    using var_sub = var_value<value_type_t<vari_sub>>;
    return var_sub(new vari_sub(vi_->bottomRows(n)));
  }

  /**
   * Return a block consisting of rows in the middle.
   * @param start_row Starting row index
   * @param n Number of rows
   */
  inline auto middleRows(Eigen::Index start_row, Eigen::Index n) const {
    using vari_sub = decltype(vi_->middleRows(start_row, n));
    using var_sub = var_value<value_type_t<vari_sub>>;
    return var_sub(new vari_sub(vi_->middleRows(start_row, n)));
  }
  inline auto middleRows(Eigen::Index start_row, Eigen::Index n) {
    using vari_sub = decltype(vi_->middleRows(start_row, n));
    using var_sub = var_value<value_type_t<vari_sub>>;
    return var_sub(new vari_sub(vi_->middleRows(start_row, n)));
  }

  /**
   * Return a block consisting of the left-most columns
   * @param n Number of columns
   */
  inline auto leftCols(Eigen::Index n) const {
    using vari_sub = decltype(vi_->leftCols(n));
    using var_sub = var_value<value_type_t<vari_sub>>;
    return var_sub(new vari_sub(vi_->leftCols(n)));
  }
  inline auto leftCols(Eigen::Index n) {
    using vari_sub = decltype(vi_->leftCols(n));
    using var_sub = var_value<value_type_t<vari_sub>>;
    return var_sub(new vari_sub(vi_->leftCols(n)));
  }

  /**
   * Return a block consisting of the right-most columns
   * @param n Number of columns
   */
  inline auto rightCols(Eigen::Index n) const {
    using vari_sub = decltype(vi_->rightCols(n));
    using var_sub = var_value<value_type_t<vari_sub>>;
    return var_sub(new vari_sub(vi_->rightCols(n)));
  }
  inline auto rightCols(Eigen::Index n) {
    using vari_sub = decltype(vi_->rightCols(n));
    using var_sub = var_value<value_type_t<vari_sub>>;
    return var_sub(new vari_sub(vi_->rightCols(n)));
  }

  /**
   * Return a block consisting of columns in the middle.
   * @param start_col Starting column index
   * @param n Number of columns
   */
  inline auto middleCols(Eigen::Index start_col, Eigen::Index n) const {
    using vari_sub = decltype(vi_->middleCols(start_col, n));
    using var_sub = var_value<value_type_t<vari_sub>>;
    return var_sub(new vari_sub(vi_->middleCols(start_col, n)));
  }
  inline auto middleCols(Eigen::Index start_col, Eigen::Index n) {
    using vari_sub = decltype(vi_->middleCols(start_col, n));
    using var_sub = var_value<value_type_t<vari_sub>>;
    return var_sub(new vari_sub(vi_->middleCols(start_col, n)));
  }

  /**
   * Return an Array.
   */
  inline auto array() const {
    using vari_sub = decltype(vi_->array());
    using var_sub = var_value<value_type_t<vari_sub>>;
    return var_sub(new vari_sub(vi_->array()));
  }
  inline auto array() {
    using vari_sub = decltype(vi_->array());
    using var_sub = var_value<value_type_t<vari_sub>>;
    return var_sub(new vari_sub(vi_->array()));
  }

  /**
   * Return an Matrix.
   */
  inline auto matrix() const {
    using vari_sub = decltype(vi_->matrix());
    using var_sub = var_value<value_type_t<vari_sub>>;
    return var_sub(new vari_sub(vi_->matrix()));
  }
  inline auto matrix() {
    using vari_sub = decltype(vi_->matrix());
    using var_sub = var_value<value_type_t<vari_sub>>;
    return var_sub(new vari_sub(vi_->matrix()));
  }

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

  /**
   * Returns number of rows. Only available if `T` is a matrix.
   * @return number of rows.
   */
  template <typename U = T,
            require_any_t<is_eigen<U>, is_matrix_cl<U>>* = nullptr>
  inline auto rows() const noexcept {
    return vi_->rows();
  }

  /**
   * Returns number of columns. Only available if `T` is a matrix.
   * @return number of columns.
   */
  template <typename U = T,
            require_any_t<is_eigen<U>, is_matrix_cl<U>>* = nullptr>
  inline auto cols() const noexcept {
    return vi_->cols();
  }

  /**
   * Assignment of another plain var value, when this also contains a plain
   * type.
   * @tparam S type of the value in the `var_value` to assing
   * @param other the value to assign
   * @return this
   */
  template <typename S, require_assignable_t<value_type, S>* = nullptr,
            require_all_plain_type_t<T, S>* = nullptr>
  inline var_value<T>& operator=(const var_value<S>& other) {
    vi_ = other.vi_;
    return *this;
  }

  /**
   * Assignment of another var value, when either this or the other one does not
   * contain a plain type.
   * @tparam S type of the value in the `var_value` to assing
   * @param other the value to assign
   * @return this
   */
  template <typename S, typename T_ = T,
            require_assignable_t<value_type, S>* = nullptr,
            require_any_not_plain_type_t<T_, S>* = nullptr>
  inline var_value<T>& operator=(const var_value<S>& other) {
    arena_t<plain_type_t<T>> prev_val = vi_->val_;
    vi_->val_ = other.val();
    // no need to change any adjoints - these are just zeros before the reverse
    // pass

    reverse_pass_callback(
        [this_vi = this->vi_, other_vi = other.vi_, prev_val]() mutable {
          this_vi->val_ = prev_val;

          // we have no way of detecting aliasing between this->vi_->adj_ and
          // other.vi_->adj_, so we must copy adjoint before reseting to zero

          // we can reuse prev_val instead of allocating a new matrix
          prev_val = this_vi->adj_;
          this_vi->adj_.setZero();
          other_vi->adj_ += prev_val;
        });
    return *this;
  }

  /**
   * No-op to match with Eigen methods which call eval
   */
  template <typename T_ = T, require_plain_type_t<T_>* = nullptr>
  inline auto& eval() noexcept {
    return *this;
  }
  template <typename T_ = T, require_plain_type_t<T_>* = nullptr>
  inline const auto& eval() const noexcept {
    return *this;
  }

  /**
   * For non-plain types evaluate to the plain type
   */
  template <typename T_ = T, require_not_plain_type_t<T_>* = nullptr>
  inline auto eval() {
    return var_value<plain_type_t<T>>(*this);
  }
  template <typename T_ = T, require_not_plain_type_t<T_>* = nullptr>
  inline auto eval() const {
    return var_value<plain_type_t<T>>(*this);
  }

  /**
   * Copy assignment operator delegates to general assignment operator
   * @param other the value to assign
   * @return this
   */
  inline var_value<T>& operator=(const var_value<T>& other) {
    return operator=<T>(other);
  }
};

// For backwards compatability the default value is double
using var = var_value<double>;

}  // namespace math

/**
 * Template specialization defining the scalar type of
 * values stored in var_value.
 *
 * @tparam T type to check.
 * @ingroup type_trait
 */
template <typename T>
struct scalar_type<T, std::enable_if_t<is_var<T>::value>> {
  using type
      = math::var_value<scalar_type_t<typename std::decay_t<T>::value_type>>;
};

}  // namespace stan
#endif
