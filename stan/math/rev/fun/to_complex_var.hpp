#ifndef STAN_MATH_REV_FUN_TO_COMPLEX_VAR_HPP
#define STAN_MATH_REV_FUN_TO_COMPLEX_VAR_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core/callback_vari.hpp>
#include <stan/math/rev/core/var.hpp>
#include <complex>

namespace stan {
namespace math {
inline var_complex to_complex_var(double x) {
  return var_complex{x};
}
inline var_complex to_complex_var(const var& x) {
  return make_callback_var(std::complex<double>(x.val()),
  [x](auto&& ret) mutable {
    x.adj() += ret.adj().real();
  });
}
inline var_complex to_complex_var(const std::complex<var>& x) {
  return make_callback_var(std::complex<double>(x.real().val(), x.imag().val()),
  // x should be trivially destructible
  [x](auto&& ret) mutable {
    x.real().adj() += ret.adj().real();
    x.imag().adj() += ret.adj().imag();
  });
}

inline var_complex to_complex_var(var real, var imag) {
  return make_callback_var(std::complex<double>(real.val(), imag.val()),
  // x should be trivially destructible
  [real, imag](auto&& ret) mutable {
    real.adj() += ret.adj().real();
    imag.adj() += ret.adj().imag();
  });
}

}
}

#endif