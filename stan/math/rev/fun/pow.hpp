#ifndef STAN_MATH_REV_FUN_POW_HPP
#define STAN_MATH_REV_FUN_POW_HPP

#include <stan/math/prim/core.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/copysign.hpp>
#include <stan/math/prim/fun/is_any_nan.hpp>
#include <stan/math/prim/fun/isnan.hpp>
#include <stan/math/prim/fun/is_nan.hpp>
#include <stan/math/prim/fun/pow.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/inv.hpp>
#include <stan/math/rev/fun/inv_sqrt.hpp>
#include <stan/math/rev/fun/inv_square.hpp>
#include <stan/math/rev/fun/is_nan.hpp>
#include <stan/math/rev/fun/log.hpp>
#include <stan/math/rev/fun/sqrt.hpp>
#include <stan/math/rev/fun/square.hpp>
#include <stan/math/rev/fun/value_of_rec.hpp>
#include <cmath>
#include <complex>
#include <type_traits>

namespace stan {
namespace math {

/**
 * Return the base raised to the power of the exponent (cmath).
 *
 * The partial derivatives are
 *
 * \f$\frac{\partial}{\partial x} \mbox{pow}(x, y) = y x^{y-1}\f$, and
 *
 * \f$\frac{\partial}{\partial y} \mbox{pow}(x, y) = x^y \ \log x\f$.
 *
 *
   \f[
   \mbox{pow}(x, y) =
   \begin{cases}
     x^y & \mbox{if } -\infty\leq x, y \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{pow}(x, y)}{\partial x} =
   \begin{cases}
     yx^{y-1} & \mbox{if } -\infty\leq x\leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{pow}(x, y)}{\partial y} =
   \begin{cases}
     x^y\ln x & \mbox{if } -\infty\leq x\leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @param base Base variable.
 * @param exponent Exponent variable.
 * @return Base raised to the exponent.
 */
inline var pow(const var& base, const var& exponent) {
  return make_callback_var(std::pow(base.val(), exponent.val()),
                           [base, exponent](auto&& vi) mutable {
                             if (base.val() == 0.0) {
                               return;  // partials zero, avoids 0 & log(0)
                             }
                             const double vi_mul = vi.adj() * vi.val();
                             base.adj() += vi_mul * exponent.val() / base.val();
                             exponent.adj() += vi_mul * std::log(base.val());
                           });
}

/**
 * Return the base raised to the power of the exponent (cmath). For matrices
 * this is performed elementwise.
 * @tparam Mat1 An Eigen type deriving from Eigen::EigenBase, a standard vector, or
 *  a `var_value` with inner Eigen type as defined above. The `scalar_type`
 *  must be a `var`.
 * @tparam Mat2 An Eigen type deriving from Eigen::EigenBase, a standard vector, or
 *  a `var_value` with inner Eigen type as defined above. The `scalar_type`
 *  must be a `var`.
 * @param base Base variable.
 * @param exponent Exponent variable.
 * @return Base raised to the exponent.
 */
template <typename Mat1, typename Mat2,
          require_all_st_var_or_arithmetic<Mat1, Mat2>* = nullptr,
          require_any_matrix_st<is_var, Mat1, Mat2>* = nullptr,
          require_all_not_stan_scalar_t<Mat1, Mat2>* = nullptr>
inline auto pow(const Mat1& base, const Mat2& exponent) {
  using val_type = decltype(
      as_array_or_scalar(value_of(base)).pow(as_array_or_scalar(value_of(exponent))).matrix().eval());
  using ret_type = return_var_matrix_t<val_type, Mat1, Mat2>;
  if (!is_constant<Mat1>::value && !is_constant<Mat2>::value) {
    using base_t = decltype(as_array_or_scalar(base));
    arena_t<promote_scalar_t<var, base_t>> arena_base = as_array_or_scalar(base);
    using exp_t = decltype(as_array_or_scalar(exponent));
    arena_t<promote_scalar_t<var, exp_t>> arena_exponent = as_array_or_scalar(exponent);
    arena_t<ret_type> ret
        = arena_base.val().pow(arena_exponent.val()).matrix();
    reverse_pass_callback([arena_base, arena_exponent, ret]() mutable {
      auto are_vals_zero = (arena_base.val() != 0.0).eval();
      auto ret_mul = (ret.adj().array() * ret.val().array()).eval();
      arena_base.adj()
          += (are_vals_zero)
                 .select(ret_mul * arena_exponent.val()
                             / arena_base.val(),
                         0);
      arena_exponent.adj()
          += (are_vals_zero)
                 .select(ret_mul * arena_base.val().log(), 0);
    });
    return ret_type(ret);
  } else if (!is_constant<Mat2>::value) {
    auto arena_base = to_arena(as_array_or_scalar(value_of(base)));
    using exp_t = decltype(as_array_or_scalar(exponent));
    arena_t<promote_scalar_t<var, exp_t>> arena_exponent = as_array_or_scalar(exponent);
    arena_t<ret_type> ret
        = arena_base.pow(arena_exponent.val()).matrix();
    reverse_pass_callback([arena_base, arena_exponent, ret]() mutable {
      arena_exponent.adj()
          += (arena_base != 0)
                 .select(ret.adj().array() * arena_base.log()
                             * ret.val().array(),
                         0);
    });
    return ret_type(ret);
  } else {
    using base_t = decltype(as_array_or_scalar(base));
    arena_t<promote_scalar_t<var, base_t>> arena_base = as_array_or_scalar(base);
    auto arena_exponent = to_arena(as_array_or_scalar(value_of(exponent)));
    arena_t<ret_type> ret
        = arena_base.val().pow(arena_exponent).matrix();
    reverse_pass_callback([arena_base, arena_exponent, ret]() mutable {
      arena_base.adj()
          += (arena_base.val() != 0)
                 .select(ret.adj().array() * arena_exponent
                             * ret.val().array() / arena_base.val(),
                         0);
    });
    return ret_type(ret);
  }
}

/**
 * Return the base raised to the power of the exponent (cmath). For matrices
 * this is performed elementwise.
 * @tparam Mat1 An Eigen type deriving from Eigen::EigenBase or
 *  a `var_value` with inner Eigen type as defined above. The `scalar_type`
 *  must be a `var` or Arithmetic.
 * @param base Base variable.
 * @param exponent Exponent variable.
 * @return Base raised to the exponent.
 */
template <typename Mat1,
          require_matrix_st<is_var_or_arithmetic, Mat1>* = nullptr>
inline auto pow(const Mat1& base, const var& exponent) {
  using ret_type = promote_scalar_t<var, plain_type_t<Mat1>>;
  if (!is_constant<Mat1>::value) {
    arena_t<ret_type> arena_base = base;
    arena_t<ret_type> ret
        = arena_base.val().array().pow(exponent.val()).matrix();
    reverse_pass_callback([arena_base, exponent, ret]() mutable {
      auto are_vals_zero = (arena_base.val().array() != 0.0).eval();
      auto ret_mul = (ret.adj().array() * ret.val().array()).eval();
      arena_base.adj().array()
          += (are_vals_zero)
                 .select(ret_mul * exponent.val() / arena_base.val().array(),
                         0);
      exponent.adj() += (are_vals_zero)
                            .select(ret_mul * arena_base.val().array().log(), 0)
                            .sum();
    });
    return ret_type(ret);
  } else {
    arena_t<promote_scalar_t<double, Mat1>> arena_base = value_of(base);
    arena_t<ret_type> ret = arena_base.array().pow(exponent.val()).matrix();
    reverse_pass_callback([arena_base, exponent, ret]() mutable {
      exponent.adj() += (arena_base.array() != 0)
                            .select(ret.adj().array() * arena_base.array().log()
                                        * ret.val().array(),
                                    0)
                            .sum();
    });
    return ret_type(ret);
  }
}

/**
 * Return the base variable raised to the power of the exponent
 * scalar (cmath).
 *
 * The derivative for the variable is
 *
 * \f$\frac{d}{dx} \mbox{pow}(x, c) = c x^{c-1}\f$.
 *
 * The template parameters are coded as they are so that arithmetic
 * types will not be promoted into the `var` slots.
 *
 * @tparam T arithmetic type
 * @param base Base variable.
 * @param exponent Exponent scalar.
 * @return Base raised to the exponent.
 */
template <typename T, require_arithmetic_t<T>* = nullptr>
inline var pow(const var& base, T exponent) {
  if (exponent == 0.5) {
    return sqrt(base);
  } else if (exponent == 1.0) {
    return base;
  } else if (exponent == 2.0) {
    return square(base);
  } else if (exponent == -2.0) {
    return inv_square(base);
  } else if (exponent == -1.0) {
    return inv(base);
  } else if (exponent == -0.5) {
    return inv_sqrt(base);
  } else {
    return make_callback_var(
        std::pow(base.val(), exponent), [base, exponent](auto&& vi) mutable {
          if (base.val() == 0.0) {
            return;  // partials zero, avoids 0 & log(0)
          }
          base.adj() += vi.adj() * exponent * vi.val() / base.val();
        });
  }
}

/**
 * Return the base matrix variable raised to the power of the exponent
 * scalar (cmath).
 *
 * The derivative for the variable is the same as the elementwise
 *
 * \f$\frac{d}{dx} \mbox{pow}(x, c) = c x^{c-1}\f$.
 *
 * @tparam Mat An Eigen type deriving from Eigen::EigenBase or
 *  a `var_value` with inner Eigen type as defined above. The `scalar_type`
 *  must be a `var`.
 * @tparam T arithmetic type
 * @param base Base matrix.
 * @param exponent Exponent scalar.
 * @return Base raised to the exponent.
 */
template <typename Mat, typename T, require_arithmetic_t<T>* = nullptr,
          require_matrix_st<is_var, Mat>* = nullptr>
inline auto pow(const Mat& base, T exponent) {
  using ret_type = plain_type_t<Mat>;
  if (exponent == 0.5) {
    return ret_type(sqrt(base));
  } else if (exponent == 1.0) {
    return ret_type(base);
  } else if (exponent == 2.0) {
    return ret_type(square(base));
  } else if (exponent == -2.0) {
    return ret_type(inv_square(base));
  } else if (exponent == -1.0) {
    return ret_type(inv(base));
  } else if (exponent == -0.5) {
    return ret_type(inv_sqrt(base));
  } else {
    arena_t<Mat> arena_base = base;
    arena_t<Mat> ret = arena_base.val().array().pow(exponent);
    reverse_pass_callback([arena_base, exponent, ret]() mutable {
      arena_base.adj().array()
          += (arena_base.val().array() != 0.0)
                 .select(ret.adj().array() * exponent * ret.val().array()
                             / arena_base.val().array(),
                         0);
    });
    return ret_type(ret);
  }
}

/**
 * Return the base scalar raised to the power of the exponent
 * variable (cmath).
 *
 * The derivative for the variable is
 *
 * \f$\frac{d}{d y} \mbox{pow}(c, y) = c^y \log c \f$.
 *
 * The template parameters are coded as they are so that arithmetic
 * types will not be promoted into the `var` slots.
 *
 * @tparam T arithmetic type
 *
 * @param base Base scalar.
 * @param exponent Exponent variable.
 * @return Base raised to the exponent.
 */
template <typename T, require_arithmetic_t<T>* = nullptr>
inline var pow(T base, const var& exponent) {
  return make_callback_var(
      std::pow(base, exponent.val()), [base, exponent](auto&& vi) mutable {
        if (base == 0.0) {
          return;  // partials zero, avoids 0 & log(0)
        }
        exponent.adj() += vi.adj() * std::log(base) * vi.val();
      });
}

/**
 * Return the base scalar raised to the power of the exponent
 * matrix elementwise.
 *
 * The derivative for the variable is
 *
 * \f$\frac{d}{d y} \mbox{pow}(c, y) = c^y \log c \f$.
 *
 *
 * @tparam T arithmetic type
 * @tparam Mat An Eigen type deriving from Eigen::EigenBase or
 *  a `var_value` with inner Eigen type as defined above. The `scalar_type`
 * must be a `var`.
 *
 * @param base Base scalar.
 * @param exponent Exponent variable.
 * @return Base raised to the exponent.
 */
template <typename T, typename Mat, require_arithmetic_t<T>* = nullptr,
          require_matrix_st<is_var, Mat>* = nullptr>
inline auto pow(T base, const Mat& exponent) {
  using ret_type = plain_type_t<Mat>;
  arena_t<ret_type> arena_exponent = exponent;
  arena_t<ret_type> ret = exponent.val().unaryExpr(
      [base](auto&& x) { return std::pow(base, x); });
  reverse_pass_callback([base, arena_exponent, ret]() mutable {
    if (base == 0.0) {
      return;  // partials zero, avoids 0 & log(0)
    }
    arena_exponent.adj().array()
        += ret.adj().array() * std::log(base) * ret.val().array();
  });
  return ret_type(ret);
}

/**
 * Return the base scalar raised to the power of the exponent
 * matrix elementwise.
 *
 * The derivative for the variable is
 *
 * \f$\frac{d}{d y} \mbox{pow}(c, y) = c^y \log c \f$.
 *
 *
 * @tparam Mat An Eigen type deriving from Eigen::EigenBase or
 *  a `var_value` with inner Eigen type as defined above. The `scalar_type`
 * must be a `var`.
 *
 * @param base Base scalar.
 * @param exponent Exponent variable.
 * @return Base raised to the exponent.
 */
template <typename Mat, require_matrix_st<is_var_or_arithmetic, Mat>* = nullptr>
inline auto pow(var base, const Mat& exponent) {
  using ret_type = promote_scalar_t<var, plain_type_t<Mat>>;
  if (!is_constant<Mat>::value) {
    arena_t<ret_type> arena_exponent = exponent;
    arena_t<ret_type> ret = arena_exponent.val().unaryExpr(
        [base_val = base.val()](auto&& x) { return std::pow(base_val, x); });
    reverse_pass_callback([base, arena_exponent, ret]() mutable {
      if (unlikely(base.val() == 0.0)) {
        return;  // partials zero, avoids 0 & log(0)
      }
      auto ret_mul = (ret.adj().array() * ret.val().array()).eval();
      base.adj() += (ret_mul * arena_exponent.val().array() / base.val()).sum();
      arena_exponent.adj().array() += ret_mul * std::log(base.val());
    });
    return ret_type(ret);
  } else {
    arena_t<promote_scalar_t<double, ret_type>> arena_exponent
        = value_of(exponent);
    arena_t<promote_scalar_t<var, ret_type>> ret = arena_exponent.unaryExpr(
        [base_val = base.val()](auto&& x) { return std::pow(base_val, x); });
    reverse_pass_callback([base, arena_exponent, ret]() mutable {
      if (unlikely(base.val() == 0.0)) {
        return;  // partials zero, avoids 0 & log(0)
      }
      base.adj() += (ret.adj().array() * arena_exponent.array()
                     * ret.val().array() / base.val())
                        .sum();
    });
    return ret_type(ret);
  }
}

// must uniquely match all pairs of { complex<var>, complex<T>, var, T }
// with at least one var and at least one complex, where T is arithmetic:
// 1) complex<var>, complex<var>
// 2) complex<var>, complex<T>
// 3) complex<var>, var
// 4) complex<var>, T
// 5) complex<T>, complex<var>
// 6) complex<T>, var
// 7) var, complex<var>
// 8) var, complex<T>
// 9) T, complex<var>

/**
 * Return the first argument raised to the power of the second argument.
 *
 * @param x first argument
 * @param y second argument
 * @return first argument to the power of the second argument
 */
inline std::complex<var> pow(const std::complex<var>& x,
                             const std::complex<var>& y) {
  return internal::complex_pow(x, y);
}

/**
 * Return the first argument raised to the power of the second argument.
 *
 * @tparam T arithmetic type
 * @param x first argument
 * @param y second argument
 * @return first argument to the power of the second argument
 */
template <typename T, typename = require_arithmetic_t<T>>
inline std::complex<var> pow(const std::complex<var>& x,
                             const std::complex<T> y) {
  return internal::complex_pow(x, y);
}

/**
 * Return the first argument raised to the power of the second argument.
 *
 * @param x first argument
 * @param y second argument
 * @return first argument to the power of the second argument
 */
inline std::complex<var> pow(const std::complex<var>& x, const var& y) {
  return internal::complex_pow(x, y);
}

/**
 * Return the first argument raised to the power of the second argument.
 *
 * @tparam T arithmetic type
 * @param x first argument
 * @param y second argument
 * @return first argument to the power of the second argument
 */
template <typename T, typename = require_arithmetic_t<T>>
inline std::complex<var> pow(const std::complex<var>& x, T y) {
  return internal::complex_pow(x, y);
}

/**
 * Return the first argument raised to the power of the second argument.
 *
 * @tparam T arithmetic type
 * @param x first argument
 * @param y second argument
 * @return first argument to the power of the second argument
 */
template <typename T, typename = require_arithmetic_t<T>>
inline std::complex<var> pow(std::complex<T> x, const std::complex<var>& y) {
  return internal::complex_pow(x, y);
}

/**
 * Return the first argument raised to the power of the second argument.
 *
 * @tparam T arithmetic type
 * @param x first argument
 * @param y second argument
 * @return first argument to the power of the second argument
 */
template <typename T, typename = require_arithmetic_t<T>>
inline std::complex<var> pow(std::complex<T> x, const var& y) {
  return internal::complex_pow(x, y);
}

/**
 * Return the first argument raised to the power of the second argument.
 *
 * @param x first argument
 * @param y second argument
 * @return first argument to the power of the second argument
 */
inline std::complex<var> pow(const var& x, const std::complex<var>& y) {
  return internal::complex_pow(x, y);
}

/**
 * Return the first argument raised to the power of the second argument.
 *
 * @tparam T arithmetic type
 * @param x first argument
 * @param y second argument
 * @return first argument to the power of the second argument
 */
template <typename T, typename = require_arithmetic_t<T>>
inline std::complex<var> pow(const var& x, std::complex<T> y) {
  return internal::complex_pow(x, y);
}

/**
 * Return the first argument raised to the power of the second argument.
 *
 * @tparam T arithmetic type
 * @param x first argument
 * @param y second argument
 * @return first argument to the power of the second argument
 */
template <typename T, typename = require_arithmetic_t<T>>
inline std::complex<var> pow(T x, const std::complex<var>& y) {
  return internal::complex_pow(x, y);
}

/**
 * Return the first argument raised to the power of the second argument.
 *
 * Note: this overload is required because gcc still provides the
 * C++99 template function `pow(complex<T>, int)`, which introduces
 * an ambiguity.
 *
 * @param x first argument
 * @param y second argument
 * @return first argument to the power of the second argument
 */
inline std::complex<var> pow(const std::complex<var>& x, int y) {
  return internal::complex_pow(x, y);
}

}  // namespace math
}  // namespace stan
#endif
