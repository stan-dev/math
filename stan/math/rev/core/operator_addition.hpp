#ifndef STAN_MATH_REV_CORE_OPERATOR_ADDITION_HPP
#define STAN_MATH_REV_CORE_OPERATOR_ADDITION_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/rev/core/var.hpp>
#include <stan/math/rev/core/vv_vari.hpp>
#include <stan/math/rev/core/vd_vari.hpp>
#include <stan/math/prim/err/is_equal.hpp>
#include <stan/math/prim/err/check_matching_dims.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/fill.hpp>
#include <stan/math/prim/fun/is_any_nan.hpp>

namespace stan {
namespace math {

namespace internal {
class add_vv_vari final : public op_vv_vari {
 public:
  add_vv_vari(vari* avi, vari* bvi)
      : op_vv_vari(avi->val_ + bvi->val_, avi, bvi) {}
  void chain() {
    if (unlikely(is_any_nan(avi_->val_, bvi_->val_))) {
      avi_->adj_ = NOT_A_NUMBER;
      bvi_->adj_ = NOT_A_NUMBER;
    } else {
      avi_->adj_ += adj_;
      bvi_->adj_ += adj_;
    }
  }
};

class add_vd_vari final : public op_vd_vari {
 public:
  add_vd_vari(vari* avi, double b) : op_vd_vari(avi->val_ + b, avi, b) {}
  void chain() {
    if (unlikely(is_any_nan(avi_->val_, bd_))) {
      avi_->adj_ = NOT_A_NUMBER;
    } else {
      avi_->adj_ += adj_;
    }
  }
};
}  // namespace internal

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

// NEWWW


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
template <typename VarMat1, typename VarMat2, require_all_rev_matrix_t<VarMat1, VarMat2>* = nullptr>
inline auto operator+(const VarMat1& a, const VarMat2& b) {
  check_matching_dims("operator+", "a", a, "b", b);
  using ret_type = decltype(a.val() + b.val());
  promote_var_matrix_t<ret_type, VarMat1, VarMat2> ret((a.val() + b.val()).eval());
    reverse_pass_callback([ret, a, b]() mutable {
      const_cast<VarMat1&>(a).adj() += ret.adj();
      const_cast<VarMat2&>(b).adj() += ret.adj();
    });
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
template <typename Arith, typename VarMat, require_st_arithmetic<Arith>* = nullptr,
 require_rev_matrix_t<VarMat>* = nullptr>
inline auto operator+(const VarMat& a, const Arith& b) {
  if (is_eigen<Arith>::value) {
    check_matching_dims("operator+", "a", a, "b", b);
  }
  using ret_inner_type = plain_type_t<decltype((a.val().array() + as_array_or_scalar(b)).matrix())>;
  using ret_type = promote_var_matrix_t<ret_inner_type, VarMat>;
  if (is_equal(b, 0.0)) {
    return ret_type(a);
  } else {
    ret_type ret(a.val().array() + as_array_or_scalar(b));
    reverse_pass_callback([ret, a]() mutable {
      const_cast<VarMat&>(a).adj() += ret.adj();
    });
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
 template <typename Arith, typename VarMat, require_st_arithmetic<Arith>* = nullptr,
  require_rev_matrix_t<VarMat>* = nullptr>
inline auto operator+(const Arith& a, const VarMat& b) {
  return b + a;
}

// NEW 222



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
template <typename Var, typename EigMat,
 require_eigen_vt<std::is_arithmetic, EigMat>* = nullptr,
 require_var_vt<std::is_arithmetic, Var>* = nullptr>
inline auto operator+(const Var& a, const EigMat& b) {
  arena_t<promote_scalar_t<var, EigMat>>  ret(a.val() + b.array());
  reverse_pass_callback([ret, a]() mutable {
    a.adj() += ret.adj().sum();
  });
  return ret;
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
template <typename EigMat, typename Var, require_var_vt<std::is_arithmetic, Var>* = nullptr,
 require_eigen_vt<std::is_arithmetic, EigMat>* = nullptr>
inline auto operator+(const EigMat& a, const Var& b) {
  return b + a;
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
template <typename Var, typename EigMat,
 require_rev_matrix_t<EigMat>* = nullptr,
 require_var_vt<std::is_arithmetic, Var>* = nullptr>
inline auto operator+(const Var& a, const EigMat& b) {
  arena_t<EigMat> ret(a.val() + b.val().array());
  reverse_pass_callback([ret, a, b]() mutable {
    a.adj() += ret.adj().sum();
    const_cast<EigMat&>(b).adj() += ret.adj();
  });
  return ret;
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
 template <typename Var, typename EigMat,
  require_rev_matrix_t<EigMat>* = nullptr,
  require_var_vt<std::is_arithmetic, Var>* = nullptr>
 inline auto operator+(const EigMat& a, const Var& b) {
   return b + a;
 }

}  // namespace math
}  // namespace stan
#endif
