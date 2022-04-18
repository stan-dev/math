#ifndef STAN_MATH_FWD_CORE_STD_NUMERIC_LIMITS_HPP
#define STAN_MATH_FWD_CORE_STD_NUMERIC_LIMITS_HPP

#include <stan/math/fwd/core/fvar.hpp>
#include <limits>

namespace std {

// inherits documentation
template <typename T>
struct numeric_limits<stan::math::fvar<T>> {
  static const bool is_specialized{true};
  
  static const int digits{numeric_limits<double>::digits};

  static const int digits10{numeric_limits<double>::digits10};

  static const bool is_signed{numeric_limits<double>::is_signed};

  static const bool is_integer{numeric_limits<double>::is_integer};

  static const bool is_exact{numeric_limits<double>::is_exact};

  static const int radix{numeric_limits<double>::radix};

  static constexpr int max_digits10{numeric_limits<T>::max_digits10};

  static const int min_exponent{numeric_limits<double>::min_exponent};

  static const int min_exponent10{numeric_limits<double>::min_exponent10};

  static const int max_exponent{numeric_limits<double>::max_exponent};

  static const int max_exponent10{numeric_limits<double>::max_exponent10};

  static const bool has_infinity{numeric_limits<double>::has_infinity};

  static const bool has_quiet_NaN{numeric_limits<double>::has_quiet_NaN};

  static const bool
  has_signaling_NaN{numeric_limits<double>::has_signaling_NaN};
  
  static const float_denorm_style
  has_denorm{numeric_limits<double>::has_denorm};

  static const bool has_denorm_loss{numeric_limits<double>::has_denorm_loss};

  static const bool is_iec559{numeric_limits<double>::is_iec559};

  static const bool is_bounded{numeric_limits<double>::is_bounded};

  static const bool is_modulo{numeric_limits<double>::is_modulo};

  static const bool traps{numeric_limits<double>::traps};

  static const bool tinyness_before{numeric_limits<double>::tinyness_before};

  static const float_round_style
  round_style{numeric_limits<double>::round_style};

  static constexpr stan::math::fvar<T> min() {
    return numeric_limits<double>::min();
  }

  static constexpr stan::math::fvar<T> max() {
    return numeric_limits<double>::max();
  }

  static constexpr stan::math::fvar<T> epsilon() {
    return numeric_limits<double>::epsilon();
  }

  static constexpr stan::math::fvar<T> round_error() {
    return numeric_limits<double>::round_error();
  }

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
};

}  // namespace std
#endif
