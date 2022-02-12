#ifndef STAN_MATH_REV_FUN_EXP_HPP
#define STAN_MATH_REV_FUN_EXP_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/cos.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/rev/fun/is_inf.hpp>
#include <stan/math/prim/fun/isinf.hpp>
#include <stan/math/prim/fun/isfinite.hpp>
#include <stan/math/rev/fun/is_nan.hpp>
#include <stan/math/rev/fun/cos.hpp>
#include <stan/math/rev/fun/sin.hpp>
#include <stan/math/rev/functor.hpp>
#include <cmath>
#include <complex>

namespace stan {
namespace math {

/**
 * Return the exponentiation (base e) of the specified complex number.
 * @param z argument
 * @return exponentiation of argument
 */
inline std::complex<var> exp(const std::complex<var>& z) {
  return internal::complex_exp(z);
}

/**
 * Return the exponentiation of the elements of x
 *
 * @tparam T type of x
 * @param x argument
 * @return elementwise exponentiation of x
 */
template <typename T, require_var_matrix_t<T>* = nullptr>
inline auto exp(const T& x) {
  return make_callback_var(
      x.val().array().exp().matrix(), [x](const auto& vi) mutable {
        x.adj() += (vi.val().array() * vi.adj().array()).matrix();
      });
}

}  // namespace math
}  // namespace stan
#endif
