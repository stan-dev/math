#ifndef STAN_MATH_FWD_SCAL_FUN_ABS_HPP
#define STAN_MATH_FWD_SCAL_FUN_ABS_HPP

#include <stan/math/fwd/cplx.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/scal/fun/value_of.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/scal/fun/abs.hpp>
#include <stan/math/prim/scal/fun/value_of.hpp>
#include <stan/math/prim/scal/meta/likely.hpp>

namespace stan {
namespace math {

template <typename T>
inline fvar<T> abs(const fvar<T>& x) {
  if (x.val_ > 0.0)
    return x;
  else if (x.val_ < 0.0)
    return fvar<T>(-x.val_, -x.d_);
  else if (x.val_ == 0.0)
    return fvar<T>(0, 0);
  else
    return fvar<T>(abs(x.val_), NOT_A_NUMBER);
}

}  // namespace math
}  // namespace stan

namespace std {
// constrained complex overload to forward zeroing fvar to std::abs
template <class T>
inline std::complex<stan::math::fvar<T>>
abs(std::complex<stan::math::fvar<T>> const& t) {
 return abs(
  std::complex<stan::math::zeroing<stan::math::fvar<T>>>(t)
 );
}
}  // namespace std
#endif
