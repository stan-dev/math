#ifndef STAN_MATH_REV_CORE_OPERATOR_SUBTRACTION_HPP
#define STAN_MATH_REV_CORE_OPERATOR_SUBTRACTION_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core/var.hpp>
#include <stan/math/rev/core/arena_matrix.hpp>
#include <stan/math/rev/core/callback_vari.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/is_any_nan.hpp>

namespace stan {
namespace math {

/**
 * Subtraction operator for variables.
 *
 * The partial derivatives are defined by
 *
 * \f$\frac{\partial}{\partial x} (x-y) = 1\f$, and
 *
 * \f$\frac{\partial}{\partial y} (x-y) = -1\f$.
 *
   \f[
   \mbox{operator-}(x, y) =
   \begin{cases}
     x-y & \mbox{if } -\infty\leq x, y \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{operator-}(x, y)}{\partial x} =
   \begin{cases}
     1 & \mbox{if } -\infty\leq x, y \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{operator-}(x, y)}{\partial y} =
   \begin{cases}
     -1 & \mbox{if } -\infty\leq x, y \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @tparam Var1 value type of a var
 * @tparam Var2 value type of a var
 * @param a First variable operand.
 * @param b Second variable operand.
 * @return Variable result of subtracting the second variable from
 * the first.
 */
inline var operator-(const var& a, const var& b) {
  return make_callback_vari(a.vi_->val_ - b.vi_->val_,
                            [avi = a.vi_, bvi = b.vi_](const auto& vi) mutable {
                              avi->adj_ += vi.adj_;
                              bvi->adj_ -= vi.adj_;
                            });
}

/**
 * Subtraction operator for variable and scalar.
 *
 * The derivative for the variable is
 *
 * \f$\frac{\partial}{\partial x} (x-c) = 1\f$, and
 *
 * @tparam Var value type of a var
 * @tparam Arith An arithmetic type
 * @param a First variable operand.
 * @param b Second scalar operand.
 * @return Result of subtracting the scalar from the variable.
 */
template <typename Arith, require_arithmetic_t<Arith>* = nullptr>
inline var operator-(const var& a, Arith b) {
  if (unlikely(b == 0.0)) {
    return a;
  }
  return make_callback_vari(
      a.vi_->val_ - b,
      [avi = a.vi_](const auto& vi) mutable { avi->adj_ += vi.adj_; });
}

/**
 * Subtraction operator for scalar and variable.
 *
 * The derivative for the variable is
 *
 * \f$\frac{\partial}{\partial y} (c-y) = -1\f$, and
 *
 * @tparam Var value type of a var
 * @tparam Arith An arithmetic type
 * @param a First scalar operand.
 * @param b Second variable operand.
 * @return Result of subtracting a variable from a scalar.
 */
template <typename Arith, require_arithmetic_t<Arith>* = nullptr>
inline var operator-(Arith a, const var& b) {
  return make_callback_vari(
      a - b.vi_->val_,
      [bvi = b.vi_, a](const auto& vi) mutable { bvi->adj_ -= vi.adj_; });
}

/**
 * Subtraction operator for matrix variables.
 *
 * @tparam VarMat1 A matrix of vars or a var with an underlying matrix type.
 * @tparam VarMat2 A matrix of vars or a var with an underlying matrix type.
 * @param a First variable operand.
 * @param b Second variable operand.
 * @return Variable result of subtracting two variables.
 */
template <typename VarMat1, typename VarMat2,
          require_all_rev_matrix_t<VarMat1, VarMat2>* = nullptr>
inline auto subtract(const VarMat1& a, const VarMat2& b) {
  check_matching_dims("subtract", "a", a, "b", b);
  using op_ret_type = decltype(a.val() - b.val());
  using ret_type = return_var_matrix_t<op_ret_type, VarMat1, VarMat2>;
  arena_t<VarMat1> arena_a = a;
  arena_t<VarMat2> arena_b = b;
  arena_t<ret_type> ret((arena_a.val() - arena_b.val()));
  reverse_pass_callback([ret, arena_a, arena_b]() mutable {
    for (Eigen::Index j = 0; j < ret.cols(); ++j) {
      for (Eigen::Index i = 0; i < ret.rows(); ++i) {
        const auto ref_adj = ret.adj().coeffRef(i, j);
        arena_a.adj().coeffRef(i, j) += ref_adj;
        arena_b.adj().coeffRef(i, j) -= ref_adj;
      }
    }
  });
  return ret_type(ret);
}

/**
 * Subtraction operator for a matrix variable and arithmetic.
 *
 * @tparam VarMat A matrix of vars or a var with an underlying matrix type.
 * @tparam Arith A type with an arithmetic Scalar type.
 * @param a First variable operand.
 * @param b Second variable operand.
 * @return Variable result of subtracting two variables.
 */
template <typename Arith, typename VarMat,
          require_st_arithmetic<Arith>* = nullptr,
          require_rev_matrix_t<VarMat>* = nullptr>
inline auto subtract(const VarMat& a, const Arith& b) {
  if (is_eigen<Arith>::value) {
    check_matching_dims("subtract", "a", a, "b", b);
  }
  using op_ret_type = plain_type_t<decltype(
      (a.val().array() - as_array_or_scalar(b)).matrix())>;
  using ret_type = return_var_matrix_t<op_ret_type, VarMat>;
  arena_t<VarMat> arena_a = a;
  arena_t<ret_type> ret(arena_a.val().array() - as_array_or_scalar(b));
  reverse_pass_callback(
      [ret, arena_a]() mutable { arena_a.adj() += ret.adj(); });
  return ret_type(ret);
}

/**
 * Subtraction operator for an arithmetic type and matrix variable.
 *
 * @tparam VarMat A matrix of vars or a var with an underlying matrix type.
 * @tparam Arith A type with an arithmetic Scalar type.
 * @param a First variable operand.
 * @param b Second variable operand.
 * @return Variable result of subtracting two variables.
 */
template <typename Arith, typename VarMat,
          require_st_arithmetic<Arith>* = nullptr,
          require_rev_matrix_t<VarMat>* = nullptr>
inline auto subtract(const Arith& a, const VarMat& b) {
  if (is_eigen<Arith>::value) {
    check_matching_dims("subtract", "a", a, "b", b);
  }
  using op_ret_type = plain_type_t<decltype(
      (as_array_or_scalar(a) - b.val().array()).matrix())>;
  using ret_type = return_var_matrix_t<op_ret_type, VarMat>;
  arena_t<VarMat> arena_b = b;
  arena_t<ret_type> ret(as_array_or_scalar(a) - arena_b.val().array());
  reverse_pass_callback(
      [ret, arena_b]() mutable { arena_b.adj() -= ret.adj_op(); });
  return ret_type(ret);
}

/**
 * Subtraction operator for an arithmetic matrix and variable.
 *
 * @tparam Var A `var_value` with an underlying arithmetic type.
 * @tparam EigMat An Eigen Matrix type with an arithmetic Scalar type.
 * @param a First variable operand.
 * @param b Second variable operand.
 * @return Variable result of subtracting two variables.
 */
template <typename Var, typename EigMat,
          require_var_vt<std::is_arithmetic, Var>* = nullptr,
          require_eigen_vt<std::is_arithmetic, EigMat>* = nullptr>
inline auto subtract(const Var& a, const EigMat& b) {
  using ret_type = return_var_matrix_t<EigMat>;
  arena_t<ret_type> ret(a.val() - b.array());
  reverse_pass_callback([ret, a]() mutable { a.adj() += ret.adj().sum(); });
  return ret_type(ret);
}

/**
 * Subtraction operator for a variable and arithmetic matrix.
 *
 * @tparam EigMat An Eigen Matrix type with an arithmetic Scalar type.
 * @tparam Var A `var_value` with an underlying arithmetic type.
 * @param a First variable operand.
 * @param b Second variable operand.
 * @return Variable result of subtracting two variables.
 */
template <typename EigMat, typename Var,
          require_eigen_vt<std::is_arithmetic, EigMat>* = nullptr,
          require_var_vt<std::is_arithmetic, Var>* = nullptr>
inline auto subtract(const EigMat& a, const Var& b) {
  using ret_type = return_var_matrix_t<EigMat>;
  arena_t<ret_type> ret(a.array() - b.val());
  reverse_pass_callback([ret, b]() mutable { b.adj() -= ret.adj().sum(); });
  return ret_type(ret);
}

/**
 * Subtraction operator for a variable and variable matrix.
 *
 * @tparam VarMat An Eigen Matrix type with a variable Scalar type or a
 * `var_value` with an underlying matrix type.
 * @tparam Var A `var_value` with an underlying arithmetic type.
 * @param a First variable operand.
 * @param b Second variable operand.
 * @return Variable result of subtracting two variables.
 */
template <typename Var, typename VarMat,
          require_var_vt<std::is_arithmetic, Var>* = nullptr,
          require_rev_matrix_t<VarMat>* = nullptr>
inline auto subtract(const Var& a, const VarMat& b) {
  using ret_type = return_var_matrix_t<VarMat>;
  arena_t<VarMat> arena_b(b);
  arena_t<ret_type> ret(a.val() - arena_b.val().array());
  reverse_pass_callback([ret, a, arena_b]() mutable {
    for (Eigen::Index j = 0; j < ret.cols(); ++j) {
      for (Eigen::Index i = 0; i < ret.rows(); ++i) {
        auto ret_adj = ret.adj().coeff(i, j);
        a.adj() += ret_adj;
        arena_b.adj().coeffRef(i, j) -= ret_adj;
      }
    }
  });
  return ret_type(ret);
}

/**
 * Subtraction operator for a variable matrix and variable.
 *
 * @tparam VarMat An Eigen Matrix type with a variable Scalar type or a
 * `var_value` with an underlying matrix type.
 * @tparam Var A `var_value` with an underlying arithmetic type.
 * @param a First variable operand.
 * @param b Second variable operand.
 * @return Variable result of subtracting two variables.
 */
template <typename Var, typename VarMat,
          require_rev_matrix_t<VarMat>* = nullptr,
          require_var_vt<std::is_arithmetic, Var>* = nullptr>
inline auto subtract(const VarMat& a, const Var& b) {
  using ret_type = return_var_matrix_t<VarMat>;
  arena_t<VarMat> arena_a(a);
  arena_t<ret_type> ret(arena_a.val().array() - b.val());
  reverse_pass_callback([ret, b, arena_a]() mutable {
    for (Eigen::Index j = 0; j < ret.cols(); ++j) {
      for (Eigen::Index i = 0; i < ret.rows(); ++i) {
        const auto ret_adj = ret.adj().coeff(i, j);
        arena_a.adj().coeffRef(i, j) += ret_adj;
        b.adj() -= ret_adj;
      }
    }
  });
  return ret_type(ret);
}

template <typename T1, typename T2,
          require_any_var_vt<std::is_arithmetic, T1, T2>* = nullptr,
          require_any_arithmetic_t<T1, T2>* = nullptr>
inline auto subtract(const T1& a, const T2& b) {
  return a - b;
}

template <typename T1, typename T2,
          require_all_var_vt<std::is_arithmetic, T1, T2>* = nullptr>
inline auto subtract(const T1& a, const T2& b) {
  return a - b;
}

/**
 * Addition operator for matrix variables.
 *
 * @tparam VarMat1 A matrix of vars or a var with an underlying matrix type.
 * @tparam VarMat2 A matrix of vars or a var with an underlying matrix type.
 * @param a First variable operand.
 * @param b Second variable operand.
 * @return Variable result of adding two variables.
 */
template <typename VarMat1, typename VarMat2,
          require_any_var_matrix_t<VarMat1, VarMat2>* = nullptr>
inline auto operator-(const VarMat1& a, const VarMat2& b) {
  return subtract(a, b);
}

}  // namespace math
}  // namespace stan
#endif
