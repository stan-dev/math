#ifndef STAN_MATH_REV_CORE_OPERATOR_DIVISION_HPP
#define STAN_MATH_REV_CORE_OPERATOR_DIVISION_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/core/operator_division.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/is_any_nan.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/rev/core/var.hpp>
#include <stan/math/rev/core/std_complex.hpp>
#include <stan/math/rev/core/operator_addition.hpp>
#include <stan/math/rev/core/operator_multiplication.hpp>
#include <stan/math/rev/core/operator_subtraction.hpp>
#include <stan/math/rev/fun/to_arena.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <complex>
#include <type_traits>

namespace stan {
namespace math {

/**
 * Division operator for two variables (C++).
 *
 * The partial derivatives for the variables are
 *
 * \f$\frac{\partial}{\partial x} (x/y) = 1/y\f$, and
 *
 * \f$\frac{\partial}{\partial y} (x/y) = -x / y^2\f$.
 *
   \f[
   \mbox{operator/}(x, y) =
   \begin{cases}
     \frac{x}{y} & \mbox{if } -\infty\leq x, y \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{operator/}(x, y)}{\partial x} =
   \begin{cases}
     \frac{1}{y} & \mbox{if } -\infty\leq x, y \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{operator/}(x, y)}{\partial y} =
   \begin{cases}
     -\frac{x}{y^2} & \mbox{if } -\infty\leq x, y \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @param dividend First variable operand.
 * @param divisor Second variable operand.
 * @return Variable result of dividing the first variable by the
 * second.
 */
inline var operator/(const var& dividend, const var& divisor) {
  return make_callback_var(
      dividend.val() / divisor.val(), [dividend, divisor](auto&& vi) {
        dividend.adj() += vi.adj() / divisor.val();
        divisor.adj()
            -= vi.adj() * dividend.val() / (divisor.val() * divisor.val());
      });
}

/**
 * Division operator for dividing a variable by a scalar (C++).
 *
 * The derivative with respect to the variable is
 *
 * \f$\frac{\partial}{\partial x} (x/c) = 1/c\f$.
 *
 * @tparam Var value type of a var
 * @tparam Arith An arithmetic type
 * @param dividend Variable operand.
 * @param divisor Scalar operand.
 * @return Variable result of dividing the variable by the scalar.
 */
template <typename Arith, require_arithmetic_t<Arith>* = nullptr>
inline var operator/(const var& dividend, Arith divisor) {
  if (divisor == 1.0) {
    return dividend;
  }
  return make_callback_var(
      dividend.val() / divisor,
      [dividend, divisor](auto&& vi) { dividend.adj() += vi.adj() / divisor; });
}

/**
 * Division operator for dividing a scalar by a variable (C++).
 *
 * The derivative with respect to the variable is
 *
 * \f$\frac{d}{d y} (c/y) = -c / y^2\f$.
 *
 * @tparam Arith An arithmetic type
 * @param dividend Scalar operand.
 * @param divisor Variable operand.
 * @return Quotient of the dividend and divisor.
 */
template <typename Arith, require_arithmetic_t<Arith>* = nullptr>
inline var operator/(Arith dividend, const var& divisor) {
  return make_callback_var(
      dividend / divisor.val(), [dividend, divisor](auto&& vi) {
        divisor.adj() -= vi.adj() * dividend / (divisor.val() * divisor.val());
      });
}

/**
 * Return matrix divided by scalar.
 *
 * @tparam Mat A type inheriting from `EigenBase` with an `Arithmetic` scalar
 * type.
 * @param[in] m specified matrix or expression
 * @param[in] c specified scalar
 * @return matrix divided by the scalar
 */
template <typename Scalar, typename Mat, require_matrix_t<Mat>* = nullptr,
          require_stan_scalar_t<Scalar>* = nullptr,
          require_all_st_var_or_arithmetic<Scalar, Mat>* = nullptr,
          require_any_st_var<Scalar, Mat>* = nullptr>
inline auto divide(const Mat& m, Scalar c) {
  if (!is_constant<Mat>::value && !is_constant<Scalar>::value) {
    arena_t<promote_scalar_t<var, Mat>> arena_m = m;
    var arena_c = c;
    auto inv_c = (1.0 / arena_c.val());
    arena_t<promote_scalar_t<var, Mat>> res = inv_c * arena_m.val();
    reverse_pass_callback([arena_c, inv_c, arena_m, res]() mutable {
      auto inv_times_adj = (inv_c * res.adj().array()).eval();
      arena_c.adj() -= (inv_times_adj * res.val().array()).sum();
      arena_m.adj().array() += inv_times_adj;
    });
    return promote_scalar_t<var, Mat>(res);
  } else if (!is_constant<Mat>::value) {
    arena_t<promote_scalar_t<var, Mat>> arena_m = m;
    auto inv_c = (1.0 / value_of(c));
    arena_t<promote_scalar_t<var, Mat>> res = inv_c * arena_m.val();
    reverse_pass_callback([inv_c, arena_m, res]() mutable {
      arena_m.adj().array() += inv_c * res.adj_op().array();
    });
    return promote_scalar_t<var, Mat>(res);
  } else {
    var arena_c = c;
    auto inv_c = (1.0 / arena_c.val());
    arena_t<promote_scalar_t<var, Mat>> res = inv_c * value_of(m).array();
    reverse_pass_callback([arena_c, inv_c, res]() mutable {
      arena_c.adj() -= inv_c * (res.adj().array() * res.val().array()).sum();
    });
    return promote_scalar_t<var, Mat>(res);
  }
}

/**
 * Return scalar divided by matrix.
 *
 * @tparam Mat Either a type inheriting from `EigenBase` with a scalar type of
 * `var` or a `var_value<T>` with type `T` inheriting from `EigenBase`.
 * @param[in] m specified matrix or expression
 * @param[in] c specified scalar
 * @return matrix divided by the scalar
 */
template <typename Scalar, typename Mat, require_matrix_t<Mat>* = nullptr,
          require_stan_scalar_t<Scalar>* = nullptr,
          require_all_st_var_or_arithmetic<Scalar, Mat>* = nullptr,
          require_any_st_var<Scalar, Mat>* = nullptr>
inline auto divide(Scalar c, const Mat& m) {
  if (!is_constant<Scalar>::value && !is_constant<Mat>::value) {
    arena_t<promote_scalar_t<var, Mat>> arena_m = m;
    auto inv_m = to_arena(arena_m.val().array().inverse());
    var arena_c = c;
    arena_t<promote_scalar_t<var, Mat>> res = arena_c.val() * inv_m;
    reverse_pass_callback([arena_c, inv_m, arena_m, res]() mutable {
      auto inv_times_res = (inv_m * res.adj().array()).eval();
      arena_m.adj().array() -= inv_times_res * res.val().array();
      arena_c.adj() += (inv_times_res).sum();
    });
    return promote_scalar_t<var, Mat>(res);
  } else if (!is_constant<Mat>::value) {
    arena_t<promote_scalar_t<var, Mat>> arena_m = m;
    auto inv_m = to_arena(arena_m.val().array().inverse());
    arena_t<promote_scalar_t<var, Mat>> res = value_of(c) * inv_m;
    reverse_pass_callback([inv_m, arena_m, res]() mutable {
      arena_m.adj().array() -= inv_m * res.adj().array() * res.val().array();
    });
    return promote_scalar_t<var, Mat>(res);
  } else {
    auto inv_m = to_arena(value_of(m).array().inverse());
    var arena_c = c;
    arena_t<promote_scalar_t<var, Mat>> res = arena_c.val() * inv_m;
    reverse_pass_callback([arena_c, inv_m, res]() mutable {
      arena_c.adj() += (inv_m * res.adj().array()).sum();
    });
    return promote_scalar_t<var, Mat>(res);
  }
}

///

/**
 * Return a matrix divided by a matrix elementwise.
 * @tparam Mat1 Either a type inheriting from `EigenBase` or a `var_value<T>`
 *  with type `T` inheriting from `EigenBase`.
 * @tparam Mat2 Either a type inheriting from `EigenBase` or a `var_value<T>`
 *  with type `T` inheriting from `EigenBase`.
 * @param[in] m1 specified matrix or expression
 * @param[in] m2 specified matrix or expression
 */
template <typename Mat1, typename Mat2,
          require_all_matrix_st<is_var_or_arithmetic, Mat1, Mat2>* = nullptr,
          require_any_matrix_st<is_var, Mat1, Mat2>* = nullptr>
inline auto divide(const Mat1& m1, const Mat2& m2) {
  if (!is_constant<Mat1>::value && !is_constant<Mat2>::value) {
    arena_t<promote_scalar_t<var, Mat1>> arena_m1 = m1;
    arena_t<promote_scalar_t<var, Mat2>> arena_m2 = m2;
    auto inv_m2 = to_arena(arena_m2.val().array().inverse());
    using val_ret = decltype((inv_m2 * arena_m1.val().array()).matrix().eval());
    using ret_type = return_var_matrix_t<val_ret, Mat1, Mat2>;
    arena_t<ret_type> res = (inv_m2.array() * arena_m1.val().array()).matrix();
    reverse_pass_callback([inv_m2, arena_m1, arena_m2, res]() mutable {
      auto inv_times_res = (inv_m2 * res.adj().array()).eval();
      arena_m1.adj().array() += inv_times_res;
      arena_m2.adj().array() -= inv_times_res * res.val().array();
    });
    return ret_type(res);
  } else if (!is_constant<Mat2>::value) {
    arena_t<promote_scalar_t<double, Mat1>> arena_m1 = value_of(m1);
    arena_t<promote_scalar_t<var, Mat2>> arena_m2 = m2;
    auto inv_m2 = to_arena(arena_m2.val().array().inverse());
    using val_ret = decltype((inv_m2 * arena_m1.array()).matrix().eval());
    using ret_type = return_var_matrix_t<val_ret, Mat1, Mat2>;
    arena_t<ret_type> res = (inv_m2.array() * arena_m1.array()).matrix();
    reverse_pass_callback([inv_m2, arena_m1, arena_m2, res]() mutable {
      arena_m2.adj().array() -= inv_m2 * res.adj().array() * res.val().array();
    });
    return ret_type(res);
  } else {
    arena_t<promote_scalar_t<var, Mat1>> arena_m1 = m1;
    arena_t<promote_scalar_t<double, Mat2>> arena_m2 = value_of(m2);
    auto inv_m2 = to_arena(arena_m2.array().inverse());
    using val_ret = decltype((inv_m2 * arena_m1.val().array()).matrix().eval());
    using ret_type = return_var_matrix_t<val_ret, Mat1, Mat2>;
    arena_t<ret_type> res = (inv_m2.array() * arena_m1.val().array()).matrix();
    reverse_pass_callback([inv_m2, arena_m1, arena_m2, res]() mutable {
      arena_m1.adj().array() += inv_m2 * res.adj().array();
    });
    return ret_type(res);
  }
}

template <typename T1, typename T2, require_any_var_matrix_t<T1, T2>* = nullptr>
inline auto operator/(const T1& dividend, const T2& divisor) {
  return divide(dividend, divisor);
}

inline std::complex<var> operator/(const std::complex<var>& x1,
                                   const std::complex<var>& x2) {
  return internal::complex_divide(x1, x2);
}

}  // namespace math
}  // namespace stan
#endif
