#ifndef STAN_MATH_FWD_CPLX_FVAR_HPP
#define STAN_MATH_FWD_CPLX_FVAR_HPP

#include <stan/math/fwd/scal/meta/ad_promotable.hpp>
#include <stan/math/fwd/scal/meta/is_fvar.hpp>
#include <stan/math/fwd/core/fvar.hpp>
#include <stan/math/prim/cplx/complex.hpp>

namespace std {

/**
 *  std::complex<fvar> inherits from stan::math::complex. Details
 *  are provided there.
 */
template <class T>
struct complex<stan::math::fvar<T>> : stan::math::complex<stan::math::fvar<T>> {
  using stan::math::complex<stan::math::fvar<T>>::complex;
};

}  // namespace std
#endif
