#ifndef STAN_MATH_REV_CORE_OPERATOR_ADDITION_HPP
#define STAN_MATH_REV_CORE_OPERATOR_ADDITION_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/prim/err/check_matching_dims.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/fill.hpp>
#include <stan/math/prim/fun/is_any_nan.hpp>
#include <stan/math/rev/core/var.hpp>
#include <stan/math/rev/core/vv_vari.hpp>
#include <stan/math/rev/core/vd_vari.hpp>

namespace stan {
namespace math {

/**
 * Addition operator for variables (C++).
 *
 * The partial derivatives are defined by
 *
 * \f$\frac{\partial}{\partial x} (x+y) = 1\f$, and
 *
 * \f$\frac{\partial}{\partial y} (x+y) = 1\f$.
 *
 *
   \f[
   \mbox{operator+}(x, y) =
   \begin{cases}
     x+y & \mbox{if } -\infty\leq x, y \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{operator+}(x, y)}{\partial x} =
   \begin{cases}
     1 & \mbox{if } -\infty\leq x, y \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{operator+}(x, y)}{\partial y} =
   \begin{cases}
     1 & \mbox{if } -\infty\leq x, y \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @param a First variable operand.
 * @param b Second variable operand.
 * @return Variable result of adding two variables.
 */
inline var operator+(const var& a, const var& b) {
  var ret(a.val() + b.val());
  if (unlikely(is_any_nan(a.val(), b.val()))) {
    reverse_pass_callback([a, b]() mutable {
      a.adj() = NOT_A_NUMBER;
      b.adj() = NOT_A_NUMBER;
    });
  } else {
    reverse_pass_callback([ret, a, b]() mutable {
      a.adj() += ret.adj();
      b.adj() += ret.adj();
    });
  }
  return ret;
}

/**
 * Addition operator for variable and scalar (C++).
 *
 * The derivative with respect to the variable is
 *
 * \f$\frac{d}{dx} (x + c) = 1\f$.
 *
 * @tparam Arith An arithmetic type
 * @param a First variable operand.
 * @param b Second scalar operand.
 * @return Result of adding variable and scalar.
 */
template <typename Arith, require_arithmetic_t<Arith>* = nullptr>
inline var operator+(const var& a, Arith b) {
  if (b == 0.0) {
    return a;
  } else {
    var ret(a.val() + b);
    if (unlikely(is_any_nan(a.val(), b))) {
      reverse_pass_callback([a]() mutable { a.adj() = NOT_A_NUMBER; });
    } else {
      reverse_pass_callback([ret, a]() mutable { a.adj() += ret.adj(); });
    }
    return ret;
  }
}

/**
 * Addition operator for scalar and variable (C++).
 *
 * The derivative with respect to the variable is
 *
 * \f$\frac{d}{dy} (c + y) = 1\f$.
 *
 * @tparam Arith An arithmetic type
 * @param a First scalar operand.
 * @param b Second variable operand.
 * @return Result of adding variable and scalar.
 */
template <typename Arith, require_arithmetic_t<Arith>* = nullptr>
inline var operator+(Arith a, const var& b) {
  return b + a;
}

/**
 * Addition operator for matrix variables (C++).
 *
 * @tparam VarMat1 A matrix of vars or a var with an underlying matrix type.
 * @tparam VarMat2 A matrix of vars or a var with an underlying matrix type.
 * @param a First variable operand.
 * @param b Second variable operand.
 * @return Variable result of adding two variables.
 */
template <typename VarMat1, typename VarMat2,
          require_all_rev_matrix_t<VarMat1, VarMat2>* = nullptr>
inline auto operator+(const VarMat1& a, const VarMat2& b) {
  check_matching_dims("operator+", "a", a, "b", b);
  using op_ret_type = decltype(a.val() + b.val());
  using ret_type = promote_var_matrix_t<op_ret_type, VarMat1, VarMat2>;
  arena_t<ret_type> ret((a.val() + b.val()).eval());
  arena_t<VarMat1> arena_a = a;
  arena_t<VarMat2> arena_b = b;
  reverse_pass_callback([ret, arena_a, arena_b]() mutable {
    for (Eigen::Index i = 0; i < ret.size(); ++i) {
      arena_a.adj().coeffRef(i) += ret.adj().coeffRef(i);
      arena_b.adj().coeffRef(i) += ret.adj().coeffRef(i);
    }
  });
  return ret_type(ret);
}

/**
 * Addition operator for a matrix variable and arithmetic (C++).
 *
 * @tparam VarMat A matrix of vars or a var with an underlying matrix type.
 * @tparam Arith A type with an arithmetic Scalar type.
 * @param a First variable operand.
 * @param b Second variable operand.
 * @return Variable result of adding two variables.
 */
template <typename Arith, typename VarMat,
          require_st_arithmetic<Arith>* = nullptr,
          require_rev_matrix_t<VarMat>* = nullptr>
inline auto operator+(const VarMat& a, const Arith& b) {
  if (is_eigen<Arith>::value) {
    check_matching_dims("operator+", "a", a, "b", b);
  }
  using op_ret_type
      = decltype((a.val().array() + as_array_or_scalar(b)).matrix());
  using ret_type = promote_var_matrix_t<op_ret_type, VarMat>;
  arena_t<VarMat> arena_a = a;
  arena_t<ret_type> ret(a.val().array() + as_array_or_scalar(b));
  reverse_pass_callback(
      [ret, arena_a]() mutable { arena_a.adj() += ret.adj_op(); });
  return ret_type(ret);
}

/**
 * Addition operator for an arithmetic type and matrix variable (C++).
 *
 * @tparam VarMat A matrix of vars or a var with an underlying matrix type.
 * @tparam Arith A type with an arithmetic Scalar type.
 * @param a First variable operand.
 * @param b Second variable operand.
 * @return Variable result of adding two variables.
 */
template <typename Arith, typename VarMat,
          require_st_arithmetic<Arith>* = nullptr,
          require_rev_matrix_t<VarMat>* = nullptr>
inline auto operator+(const Arith& a, const VarMat& b) {
  return b + a;
}

/**
 * Addition operator for an arithmetic matrix and variable (C++).
 *
 * @tparam Var A `var_value` with an underlying arithmetic type.
 * @tparam EigMat An Eigen Matrix type with an arithmetic Scalar type.
 * @param a First variable operand.
 * @param b Second variable operand.
 * @return Variable result of adding two variables.
 */
template <typename Var, typename EigMat,
          require_eigen_vt<std::is_arithmetic, EigMat>* = nullptr,
          require_var_vt<std::is_arithmetic, Var>* = nullptr>
inline auto operator+(const Var& a, const EigMat& b) {
  using ret_type = promote_scalar_t<var, EigMat>;
  arena_t<ret_type> ret(a.val() + b.array());
  reverse_pass_callback([ret, a]() mutable { a.adj() += ret.adj().sum(); });
  return ret_type(ret);
}

/**
 * Addition operator for a variable and arithmetic matrix (C++).
 *
 * @tparam EigMat An Eigen Matrix type with an arithmetic Scalar type.
 * @tparam Var A `var_value` with an underlying arithmetic type.
 * @param a First variable operand.
 * @param b Second variable operand.
 * @return Variable result of adding two variables.
 */
template <typename EigMat, typename Var,
          require_var_vt<std::is_arithmetic, Var>* = nullptr,
          require_eigen_vt<std::is_arithmetic, EigMat>* = nullptr>
inline auto operator+(const EigMat& a, const Var& b) {
  return b + a;
}

/**
 * Addition operator for a variable and variable matrix (C++).
 *
 * @tparam VarMat An Eigen Matrix type with a variable Scalar type or a
 * `var_value` with an underlying matrix type.
 * @tparam Var A `var_value` with an underlying arithmetic type.
 * @param a First variable operand.
 * @param b Second variable operand.
 * @return Variable result of adding two variables.
 */
template <typename Var, typename VarMat,
          require_rev_matrix_t<VarMat>* = nullptr,
          require_var_vt<std::is_arithmetic, Var>* = nullptr>
inline auto operator+(const Var& a, const VarMat& b) {
  arena_t<VarMat> arena_b(b);
  arena_t<VarMat> ret(a.val() + b.val().array());
  reverse_pass_callback([ret, a, arena_b]() mutable {
    for (Eigen::Index i = 0; i < arena_b.size(); ++i) {
      a.adj() += ret.adj().coeffRef(i);
      arena_b.adj().coeffRef(i) += ret.adj().coeffRef(i);
    }
  });
  return plain_type_t<VarMat>(ret);
}

/**
 * Addition operator for a variable matrix and variable (C++).
 *
 * @tparam VarMat An Eigen Matrix type with a variable Scalar type or a
 * `var_value` with an underlying matrix type.
 * @tparam Var A `var_value` with an underlying arithmetic type.
 * @param a First variable operand.
 * @param b Second variable operand.
 * @return Variable result of adding two variables.
 */
template <typename Var, typename VarMat,
          require_rev_matrix_t<VarMat>* = nullptr,
          require_var_vt<std::is_arithmetic, Var>* = nullptr>
inline auto operator+(const VarMat& a, const Var& b) {
  return b + a;
}

}  // namespace math
}  // namespace stan
#endif
