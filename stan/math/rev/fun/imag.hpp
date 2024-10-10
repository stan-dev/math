#ifndef STAN_MATH_REV_FUN_IMAG_HPP
#define STAN_MATH_REV_FUN_IMAG_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core/var.hpp>
#include <stan/math/prim/fun/imag.hpp>

namespace stan {
namespace math {
inline var imag(const var_value<std::complex<double>>& x) {
  return make_callback_var(x.val().imag(),
                           [x](auto&& ret) mutable {
                             stan::math::imag(x.adj()) += ret.adj();
                           });
}
}
}
#endif
