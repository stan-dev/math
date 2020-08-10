#ifndef STAN_MATH_REV_CORE_STD_NUMERIC_LIMITS_HPP
#define STAN_MATH_REV_CORE_STD_NUMERIC_LIMITS_HPP

#include <stan/math/rev/core/var.hpp>
#include <limits>

namespace std {

/**
 * Specialization of numeric limits for var objects.
 *
 * This implementation of std::numeric_limits<stan::math::var>
 * is used to treat var objects like value_types.
 */
template <typename T>
struct numeric_limits<stan::math::var_value<T>> {
  typedef stan::promote_args_t<T> value_type;
  static constexpr bool is_specialized = true;
  static constexpr stan::math::var_value<T> min() noexcept {
    return numeric_limits<value_type>::min();
  }
  static constexpr stan::math::var_value<T> max() noexcept {
    return numeric_limits<value_type>::max();
  }
  static constexpr int digits = numeric_limits<value_type>::digits;
  static constexpr int digits10 = numeric_limits<value_type>::digits10;
  static constexpr int max_digits10 = numeric_limits<value_type>::max_digits10;
  static constexpr bool is_signed = numeric_limits<value_type>::is_signed;
  static constexpr bool is_integer = numeric_limits<value_type>::is_integer;
  static constexpr bool is_exact = numeric_limits<value_type>::is_exact;
  static constexpr int radix = numeric_limits<value_type>::radix;
  static constexpr stan::math::var_value<T> epsilon() noexcept {
    return numeric_limits<value_type>::epsilon();
  }
  static constexpr stan::math::var_value<T> round_error() noexcept {
    return numeric_limits<value_type>::round_error();
  }
  static constexpr T lowest() noexcept {
    return numeric_limits<value_type>::lowest();
  };

  static constexpr int min_exponent = numeric_limits<value_type>::min_exponent;
  static constexpr int min_exponent10
      = numeric_limits<value_type>::min_exponent10;
  static constexpr int max_exponent = numeric_limits<value_type>::max_exponent;
  static constexpr int max_exponent10
      = numeric_limits<value_type>::max_exponent10;

  static constexpr bool has_infinity = numeric_limits<value_type>::has_infinity;
  static constexpr bool has_quiet_NaN
      = numeric_limits<value_type>::has_quiet_NaN;
  static constexpr bool has_signaling_NaN
      = numeric_limits<value_type>::has_signaling_NaN;
  static constexpr float_denorm_style has_denorm
      = numeric_limits<value_type>::has_denorm;
  static constexpr bool has_denorm_loss
      = numeric_limits<value_type>::has_denorm_loss;
  static constexpr stan::math::var_value<T> infinity() noexcept {
    return numeric_limits<value_type>::infinity();
  }
  static constexpr stan::math::var_value<T> quiet_NaN() noexcept {
    return numeric_limits<value_type>::quiet_NaN();
  }
  static constexpr stan::math::var_value<T> signaling_NaN() noexcept {
    return numeric_limits<value_type>::signaling_NaN();
  }
  static constexpr stan::math::var_value<T> denorm_min() noexcept {
    return numeric_limits<value_type>::denorm_min();
  }

  static constexpr bool is_iec559 = numeric_limits<value_type>::is_iec559;
  static constexpr bool is_bounded = numeric_limits<value_type>::is_bounded;
  static constexpr bool is_modulo = numeric_limits<value_type>::is_modulo;

  static constexpr bool traps = numeric_limits<value_type>::traps;
  static constexpr bool tinyness_before
      = numeric_limits<value_type>::tinyness_before;
  static constexpr float_round_style round_style
      = numeric_limits<value_type>::round_style;
};

}  // namespace std
#endif
