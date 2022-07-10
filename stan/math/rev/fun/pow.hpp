#ifndef STAN_MATH_REV_FUN_POW_HPP
#define STAN_MATH_REV_FUN_POW_HPP

#include <stan/math/prim/core.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
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
template <typename Scal1, typename Scal2,
          require_any_st_var<Scal1, Scal2>* = nullptr,
          require_all_stan_scalar_t<Scal1, Scal2>* = nullptr>
inline var pow(const Scal1& base, const Scal2& exponent) {
  if (is_constant<Scal2>::value) {
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
    }
  }
  return make_callback_var(
      std::pow(value_of(base), value_of(exponent)),
      [base, exponent](auto&& vi) mutable {
        if (value_of(base) == 0.0) {
          return;  // partials zero, avoids 0 & log(0)
        }
        const double vi_mul = vi.adj() * vi.val();

        if (!is_constant<Scal1>::value) {
          forward_as<var>(base).adj()
              += vi_mul * value_of(exponent) / value_of(base);
        }
        if (!is_constant<Scal2>::value) {
          forward_as<var>(exponent).adj() += vi_mul * std::log(value_of(base));
        }
      });
}

/**
 * Return the base raised to the power of the exponent (cmath). For matrices
 * this is performed elementwise.
 * @tparam Mat1 An Eigen type deriving from Eigen::EigenBase, a standard vector,
 * or a `var_value` with inner Eigen type as defined above. The `scalar_type`
 *  must be a `var`.
 * @tparam Mat2 An Eigen type deriving from Eigen::EigenBase, a standard vector,
 * or a `var_value` with inner Eigen type as defined above. The `scalar_type`
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
  check_consistent_sizes("pow", "base", base, "exponent", exponent);

  using val_type = decltype(as_array_or_scalar(value_of(base))
                                .pow(as_array_or_scalar(value_of(exponent)))
                                .matrix()
                                .eval());
  using ret_type = return_var_matrix_t<val_type, Mat1, Mat2>;
  using base_t = decltype(as_array_or_scalar(base));
  using exp_t = decltype(as_array_or_scalar(exponent));
  using base_arena_t = arena_t<base_t>;
  using exp_arena_t = arena_t<exp_t>;

  base_arena_t arena_base = as_array_or_scalar(base);
  exp_arena_t arena_exponent = as_array_or_scalar(exponent);
  arena_t<ret_type> ret
      = value_of(arena_base).pow(value_of(arena_exponent)).matrix();

  reverse_pass_callback([arena_base, arena_exponent, ret]() mutable {
    const auto& are_vals_zero = to_ref(value_of(arena_base) != 0.0);
    const auto& ret_mul = to_ref(ret.adj().array() * ret.val().array());
    if (!is_constant<Mat1>::value) {
      using base_var_arena_t = arena_t<promote_scalar_t<var, base_arena_t>>;
      forward_as<base_var_arena_t>(arena_base).adj()
          += (are_vals_zero)
                 .select(
                     ret_mul * value_of(arena_exponent) / value_of(arena_base),
                     0);
    }
    if (!is_constant<Mat2>::value) {
      using exp_var_arena_t = arena_t<promote_scalar_t<var, exp_arena_t>>;
      forward_as<exp_var_arena_t>(arena_exponent).adj()
          += (are_vals_zero).select(ret_mul * value_of(arena_base).log(), 0);
    }
  });
  return ret_type(ret);
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
template <typename Mat1, typename Scal1,
          require_all_st_var_or_arithmetic<Mat1, Scal1>* = nullptr,
          require_all_matrix_st<is_var, Mat1>* = nullptr,
          require_stan_scalar_t<Scal1>* = nullptr>
inline auto pow(const Mat1& base, const Scal1& exponent) {
  using ret_type = promote_scalar_t<var, plain_type_t<Mat1>>;

  if (is_constant<Scal1>::value) {
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
    }
  }

  arena_t<plain_type_t<Mat1>> arena_base = base;
  arena_t<ret_type> ret
      = value_of(arena_base).array().pow(value_of(exponent)).matrix();

  reverse_pass_callback([arena_base, exponent, ret]() mutable {
    const auto& are_vals_zero = to_ref(value_of(arena_base).array() != 0.0);
    const auto& ret_mul = to_ref(ret.adj().array() * ret.val().array());
    if (!is_constant<Mat1>::value) {
      forward_as<ret_type>(arena_base).adj().array()
          += (are_vals_zero)
                 .select(ret_mul * value_of(exponent)
                             / value_of(arena_base).array(),
                         0);
    }
    if (!is_constant<Scal1>::value) {
      forward_as<var>(exponent).adj()
          += (are_vals_zero)
                 .select(ret_mul * value_of(arena_base).array().log(), 0)
                 .sum();
    }
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
template <typename Scal1, typename Mat1,
          require_all_st_var_or_arithmetic<Scal1, Mat1>* = nullptr,
          require_stan_scalar_t<Scal1>* = nullptr,
          require_all_matrix_st<is_var, Mat1>* = nullptr>
inline auto pow(Scal1 base, const Mat1& exponent) {
  using ret_type = promote_scalar_t<var, plain_type_t<Mat1>>;
  arena_t<Mat1> arena_exponent = exponent;
  arena_t<ret_type> ret
      = Eigen::pow(value_of(base), value_of(arena_exponent).array());

  reverse_pass_callback([base, arena_exponent, ret]() mutable {
    if (unlikely(value_of(base) == 0.0)) {
      return;  // partials zero, avoids 0 & log(0)
    }
    const auto& ret_mul = to_ref(ret.adj().array() * ret.val().array());
    if (!is_constant<Scal1>::value) {
      forward_as<var>(base).adj()
          += (ret_mul * value_of(arena_exponent).array() / value_of(base))
                 .sum();
    }
    if (!is_constant<Mat1>::value) {
      forward_as<ret_type>(arena_exponent).adj().array()
          += ret_mul * std::log(value_of(base));
    }
  });
  return ret_type(ret);
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
