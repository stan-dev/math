#ifndef STAN_MATH_REV_FUN_FREXP_HPP
#define STAN_MATH_REV_FUN_FREXP_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>

namespace stan {
namespace math {
inline auto frexp(stan::math::var x, int* exponent) noexcept {
  return std::frexp(x.val(), exponent);
}
}
}
#endif
