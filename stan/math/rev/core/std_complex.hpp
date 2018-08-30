#ifndef STAN_MATH_REV_CORE_STD_COMPLEX_HPP
#define STAN_MATH_REV_CORE_STD_COMPLEX_HPP

#include <stan/math/rev/core/var.hpp>
#include <stan/math/rev/core/std_numeric_limits.hpp>
#include <complex>

namespace std {

/**
 * Specialization of complex for var objects.
 *
 */
template <>
constexpr complex<stan::math::var>::complex(const stan::math::var& re,
                                            const stan::math::var& im) {
  real(stan::math::is_uninitialized(re) ? stan::math::var{0.0} : re);
  imag(stan::math::is_uninitialized(im) ? stan::math::var{0.0} : im);
}

// struct complex<stan::math::var> {
//   complex()
//       : _M_real(0), _M_imag(0) {
//   }
// };

}  // namespace std
#endif
