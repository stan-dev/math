#ifndef STAN_MATH_REV_CORE_OPERATOR_SUBTRACTION_HPP
#define STAN_MATH_REV_CORE_OPERATOR_SUBTRACTION_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/is_equal.hpp>
#include <stan/math/prim/err/check_matching_dims.hpp>
#include <stan/math/rev/core/var.hpp>
#include <stan/math/rev/core/vv_vari.hpp>
#include <stan/math/rev/core/vd_vari.hpp>
#include <stan/math/rev/core/dv_vari.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/is_any_nan.hpp>

namespace stan {
namespace math {

namespace internal {
class subtract_vv_vari final : public op_vv_vari {
 public:
  subtract_vv_vari(vari* avi, vari* bvi)
      : op_vv_vari(avi->val_ - bvi->val_, avi, bvi) {}
  void chain() {
    if (unlikely(is_any_nan(avi_->val_, bvi_->val_))) {
      avi_->adj_ = NOT_A_NUMBER;
      bvi_->adj_ = NOT_A_NUMBER;
    } else {
      avi_->adj_ += adj_;
      bvi_->adj_ -= adj_;
    }
  }
};

class subtract_vd_vari final : public op_vd_vari {
 public:
  subtract_vd_vari(vari* avi, double b) : op_vd_vari(avi->val_ - b, avi, b) {}
  void chain() {
    if (unlikely(is_any_nan(avi_->val_, bd_))) {
      avi_->adj_ = NOT_A_NUMBER;
    } else {
      avi_->adj_ += adj_;
    }
  }
};

class subtract_dv_vari final : public op_dv_vari {
 public:
  subtract_dv_vari(double a, vari* bvi) : op_dv_vari(a - bvi->val_, a, bvi) {}
  void chain() {
    if (unlikely(is_any_nan(ad_, bvi_->val_))) {
      bvi_->adj_ = NOT_A_NUMBER;
    } else {
      bvi_->adj_ -= adj_;
    }
  }
};
}  // namespace internal

/**
 * Subtraction operator for variables (C++).
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
   var ret(a.val() - b.val());
   if (unlikely(is_any_nan(a.val(), b.val()))) {
     reverse_pass_callback([a, b]() mutable {
         a.adj() = NOT_A_NUMBER;
         b.adj() = NOT_A_NUMBER;
       });
   } else {
     reverse_pass_callback([ret, a, b]() mutable {
       a.adj() += ret.adj();
       b.adj() -= ret.adj();
     });
   }
   return ret;
 }

 /**
  * Subtraction operator for variable and scalar (C++).
  *
  * The derivative with respect to the variable is
  *
  * \f$\frac{d}{dx} (x + c) = 1\f$.
  *
  * @tparam Arith An arithmetic type
  * @param a First variable operand.
  * @param b Second scalar operand.
  * @return Result of subtracting variable and scalar.
  */
 template <typename Arith, require_arithmetic_t<Arith>* = nullptr>
 inline var operator-(const var& a, Arith b) {
   if (b == 0.0) {
     return a;
   } else {
     var ret(a.val() - b);
     if (unlikely(is_any_nan(a.val(), b))) {
       reverse_pass_callback([a]() mutable {
           a.adj() = NOT_A_NUMBER;
         });
     } else {
       reverse_pass_callback([ret, a]() mutable {
         a.adj() += ret.adj();
       });
     }
     return ret;
   }
 }

 /**
  * Subtraction operator for scalar and variable (C++).
  *
  * The derivative with respect to the variable is
  *
  * \f$\frac{d}{dy} (c + y) = 1\f$.
  *
  * @tparam Arith An arithmetic type
  * @param a First scalar operand.
  * @param b Second variable operand.
  * @return Result of subtracting variable and scalar.
  */
 template <typename Arith, require_arithmetic_t<Arith>* = nullptr>
 inline var operator-(Arith a, const var& b) {
     var ret(a - b.val());
     if (unlikely(is_any_nan(a, b.val()))) {
       reverse_pass_callback([b]() mutable {
           b.adj() = NOT_A_NUMBER;
         });
     } else {
       reverse_pass_callback([ret, b]() mutable {
         b.adj() -= ret.adj();
       });
     }
     return ret;
 }

 /**
  * Subtraction operator for matrix variables (C++).
  *
  * @tparam VarMat1 A matrix of vars or a var with an underlying matrix type.
  * @tparam VarMat2 A matrix of vars or a var with an underlying matrix type.
  * @param a First variable operand.
  * @param b Second variable operand.
  * @return Variable result of subtracting two variables.
  */
 template <typename VarMat1, typename VarMat2, require_all_rev_matrix_t<VarMat1, VarMat2>* = nullptr>
 inline auto operator-(const VarMat1& a, const VarMat2& b) {
   check_matching_dims("operator-", "a", a, "b", b);
   using ret_type = decltype(a.val() - b.val());
   promote_var_matrix_t<ret_type, VarMat1, VarMat2> ret((a.val() - b.val()).eval());
   arena_t<VarMat1> arena_a = a;
   arena_t<VarMat2> arena_b = b;
   reverse_pass_callback([ret, arena_a, arena_b]() mutable {
       arena_a.adj() += ret.adj_op();
       arena_b.adj() -= ret.adj_op();
     });
   return ret;
 }

 /**
  * Subtraction operator for a matrix variable and arithmetic (C++).
  *
  * @tparam VarMat A matrix of vars or a var with an underlying matrix type.
  * @tparam Arith A type with an arithmetic Scalar type.
  * @param a First variable operand.
  * @param b Second variable operand.
  * @return Variable result of subtracting two variables.
  */
 template <typename Arith, typename VarMat, require_st_arithmetic<Arith>* = nullptr,
  require_rev_matrix_t<VarMat>* = nullptr>
 inline auto operator-(const VarMat& a, const Arith& b) {
   if (is_eigen<Arith>::value) {
     check_matching_dims("operator-", "a", a, "b", b);
   }
   using ret_inner_type = plain_type_t<decltype((a.val().array() - as_array_or_scalar(b)).matrix())>;
   using ret_type = promote_var_matrix_t<ret_inner_type, VarMat>;
   if (is_equal(b, 0.0)) {
     return ret_type(a);
   } else {
     arena_t<VarMat> arena_a = a;
     ret_type ret(a.val().array() - as_array_or_scalar(b));
     reverse_pass_callback([ret, arena_a]() mutable {
       arena_a.adj() += ret.adj_op();
     });
     return ret;
   }
 }

 /**
  * Subtraction operator for an arithmetic type and matrix variable (C++).
  *
  * @tparam VarMat A matrix of vars or a var with an underlying matrix type.
  * @tparam Arith A type with an arithmetic Scalar type.
  * @param a First variable operand.
  * @param b Second variable operand.
  * @return Variable result of subtracting two variables.
  */
  template <typename Arith, typename VarMat, require_st_arithmetic<Arith>* = nullptr,
   require_rev_matrix_t<VarMat>* = nullptr>
 inline auto operator-(const Arith& a, const VarMat& b) {
   if (is_eigen<Arith>::value) {
     check_matching_dims("operator-", "a", a, "b", b);
   }
   using ret_inner_type = plain_type_t<decltype((as_array_or_scalar(a) - b.val().array()).matrix())>;
   using ret_type = promote_var_matrix_t<ret_inner_type, VarMat>;
   arena_t<VarMat> arena_b = b;
   ret_type ret(as_array_or_scalar(a) - b.val().array());
   reverse_pass_callback([ret, arena_b]() mutable {
     arena_b.adj() -= ret.adj_op();
   });
   return ret;
 }

 /**
  * Subtraction operator for an arithmetic matrix and variable (C++).
  *
  * @tparam Var A `var_value` with an underlying arithmetic type.
  * @tparam EigMat An Eigen Matrix type with an arithmetic Scalar type.
  * @param a First variable operand.
  * @param b Second variable operand.
  * @return Variable result of subtracting two variables.
  */
 template <typename Var, typename EigMat,
  require_eigen_vt<std::is_arithmetic, EigMat>* = nullptr,
  require_var_vt<std::is_arithmetic, Var>* = nullptr>
 inline auto operator-(const Var& a, const EigMat& b) {
   arena_t<promote_scalar_t<var, EigMat>>  ret(a.val() - b.array());
   reverse_pass_callback([ret, a]() mutable {
     a.adj() += ret.adj().sum();
   });
   return ret;
 }


 /**
  * Subtraction operator for a variable and arithmetic matrix (C++).
  *
  * @tparam EigMat An Eigen Matrix type with an arithmetic Scalar type.
  * @tparam Var A `var_value` with an underlying arithmetic type.
  * @param a First variable operand.
  * @param b Second variable operand.
  * @return Variable result of subtracting two variables.
  */
 template <typename EigMat, typename Var, require_var_vt<std::is_arithmetic, Var>* = nullptr,
  require_eigen_vt<std::is_arithmetic, EigMat>* = nullptr>
 inline auto operator-(const EigMat& a, const Var& b) {
   arena_t<promote_scalar_t<var, EigMat>>  ret(a.array() - b.val());
   reverse_pass_callback([ret, b]() mutable {
     b.adj() -= ret.adj().sum();
   });
   return ret;
 }

 /**
  * Subtraction operator for a variable and variable matrix (C++).
  *
  * @tparam VarMat An Eigen Matrix type with a variable Scalar type or a `var_value` with an underlying matrix type.
  * @tparam Var A `var_value` with an underlying arithmetic type.
  * @param a First variable operand.
  * @param b Second variable operand.
  * @return Variable result of subtracting two variables.
  */
 template <typename Var, typename VarMat,
  require_rev_matrix_t<VarMat>* = nullptr,
  require_var_vt<std::is_arithmetic, Var>* = nullptr>
 inline auto operator-(const Var& a, const VarMat& b) {
   arena_t<VarMat> arena_b(b);
   arena_t<VarMat> ret(a.val() - b.val().array());
   reverse_pass_callback([ret, a, arena_b]() mutable {
     a.adj() += ret.adj().sum();
     arena_b.adj() -= ret.adj();
   });
   return ret;
 }


 /**
  * Subtraction operator for a variable matrix and variable (C++).
  *
  * @tparam VarMat An Eigen Matrix type with a variable Scalar type or a `var_value` with an underlying matrix type.
  * @tparam Var A `var_value` with an underlying arithmetic type.
  * @param a First variable operand.
  * @param b Second variable operand.
  * @return Variable result of subtracting two variables.
  */
  template <typename Var, typename VarMat,
   require_rev_matrix_t<VarMat>* = nullptr,
   require_var_vt<std::is_arithmetic, Var>* = nullptr>
  inline auto operator-(const VarMat& a, const Var& b) {
    arena_t<VarMat> arena_a(a);
    arena_t<VarMat> ret(a.val().array() - b.val());
    reverse_pass_callback([ret, b, arena_a]() mutable {
      arena_a.adj().array() += ret.adj().array();
      b.adj() -= ret.adj().sum();
    });
    return ret;
  }

 }  // namespace math
 }  // namespace stan
 #endif
