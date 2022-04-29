#ifndef STAN_MATH_REV_FUN_SIGN_HPP
#define STAN_MATH_REV_FUN_SIGN_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>

namespace stan {
namespace math {
inline int sign(stan::math::var z) {
  return (z == 0) ? 0 : z < 0 ? -1 : 1;
}
}
}
#endif
