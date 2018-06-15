#ifndef STAN_MATH_REV_CPLX_VAR_HPP
#define STAN_MATH_REV_CPLX_VAR_HPP

#include <stan/math/rev/scal/meta/is_var.hpp>
#include <stan/math/rev/core/var.hpp>
#include <stan/math/prim/cplx/complex.hpp>

namespace std {

/**
 *  std::complex<var> inherits from stan::math::complex. Details
 *  are provided there.
 */
template <>
struct std::complex<stan::math::var>
    : stan::math::complex<stan::math::var> {
  using stan::math::complex<stan::math::var>::complex;
};

}  // namespace std

#endif
