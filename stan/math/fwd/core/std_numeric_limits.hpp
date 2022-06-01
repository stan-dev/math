#ifndef STAN_MATH_FWD_CORE_STD_NUMERIC_LIMITS_HPP
#define STAN_MATH_FWD_CORE_STD_NUMERIC_LIMITS_HPP

#include <stan/math/fwd/core/fvar.hpp>
#include <limits>

namespace std {

template <typename T>
struct numeric_limits<stan::math::fvar<T>> {
  static constexpr bool is_specialized{true};
  static constexpr stan::math::fvar<T> min() {
    return numeric_limits<double>::min();
  }
  static constexpr stan::math::fvar<T> max() {
    return numeric_limits<double>::max();
  }
  static constexpr int digits{numeric_limits<T>::digits};
  static constexpr int digits10{numeric_limits<T>::digits10};
  static constexpr bool is_signed{numeric_limits<T>::is_signed};
  static constexpr bool is_integer{numeric_limits<T>::is_integer};
  static constexpr bool is_exact{numeric_limits<T>::is_exact};
  static constexpr int radix{numeric_limits<T>::radix};
  static constexpr stan::math::fvar<T> epsilon() {
    return numeric_limits<double>::epsilon();
  }
  static constexpr stan::math::fvar<T> round_error() {
    return numeric_limits<double>::round_error();
  }
  static constexpr int max_digits10{numeric_limits<T>::max_digits10};
  static constexpr int min_exponent{numeric_limits<T>::min_exponent};
  static constexpr int min_exponent10{numeric_limits<T>::min_exponent10};
  static constexpr int max_exponent{numeric_limits<T>::max_exponent};
  static constexpr int max_exponent10{numeric_limits<T>::max_exponent10};

  static constexpr bool has_infinity{numeric_limits<T>::has_infinity};
  static constexpr bool has_quiet_NaN{numeric_limits<T>::has_quiet_NaN};
  static constexpr bool has_signaling_NaN{numeric_limits<T>::has_signaling_NaN};

  static constexpr float_denorm_style has_denorm{numeric_limits<T>::has_denorm};
  static constexpr bool has_denorm_loss{numeric_limits<T>::has_denorm_loss};
  static constexpr stan::math::fvar<T> infinity() {
    return numeric_limits<double>::infinity();
  }
  static constexpr stan::math::fvar<T> quiet_NaN() {
    return numeric_limits<double>::quiet_NaN();
  }
  static constexpr stan::math::fvar<T> signaling_NaN() {
    return numeric_limits<double>::signaling_NaN();
  }
  static constexpr stan::math::fvar<T> denorm_min() {
    return numeric_limits<double>::denorm_min();
  }

  static constexpr bool is_iec559{numeric_limits<T>::is_iec559};
  static constexpr bool is_bounded{numeric_limits<T>::is_bounded};
  static constexpr bool is_modulo{numeric_limits<T>::is_modulo};

  static constexpr bool traps{numeric_limits<T>::traps};
  static constexpr bool tinyness_before{numeric_limits<T>::tinyness_before};
  static constexpr float_round_style round_style{
      numeric_limits<T>::round_style};
};

}  // namespace std
#endif
