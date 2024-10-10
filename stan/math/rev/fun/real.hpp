#ifndef STAN_MATH_REV_FUN_REAL_HPP
#define STAN_MATH_REV_FUN_REAL_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core/var.hpp>
#include <stan/math/prim/fun/real.hpp>

namespace stan {
namespace math {
inline var real(const var_value<std::complex<double>>& x) {
  return make_callback_var(x.val().real(),
                           [x](auto&& ret) mutable {
                             stan::math::real(x.adj()) += ret.adj();
                           });
}
}
}
#endif
