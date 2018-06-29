#ifndef STAN_MATH_PRIM_CPLX_ZEROING_HPP
#define STAN_MATH_PRIM_CPLX_ZEROING_HPP

#include <stan/math/prim/scal/meta/is_arith_like.hpp>

namespace stan {
namespace math {

/**
 *  zeroing<T> is a zeroing T that std::complex<T> implicitly
 *  uses  when inheriting from stan::math::complex. It has the
 *  property that zeroing() results in T(0). See stan's
 *  complex class for details. Without this, in libstdc++'s
 *  complex class, when T() is called, vars (and fvar<vars>)
 *  will have a null pointer, which is bad because libstdc++
 *  immediately uses such expressions as if they were
 *  initialized to zero..
 */
template <class T = double>
struct zeroing : T {
  using T::T;  // inherit ctors
  template <std::enable_if_t<stan::is_arith_like<T>::value>* = nullptr>
  zeroing(T const& t = 0) : T(t) {}  // NOLINT(runtime/explicit)
};

}  // namespace math
}  // namespace stan
#endif
