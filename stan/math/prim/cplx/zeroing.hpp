#ifndef STAN_MATH_PRIM_CPLX_ZEROING_HPP
#define STAN_MATH_PRIM_CPLX_ZEROING_HPP

#include <stan/math/prim/scal/meta/is_arith_like.hpp>

namespace stan {
namespace math {

/**
 *  z<T> is a zeroing T that std::complex<T> implicitly uses
 *  when inheriting from stan::math::complex. It has the 
 *  property that z() results in T(0). See stan's
 *  complex class for details.
 */
template <class T = double>
struct zeroing : T {
  using T::T;  // inherit ctors
  template <std::enable_if_t<
   stan::is_arith_like<T>::value>* = nullptr>
  zeroing(T const& t =0) : T(t) {}  // NOLINT(runtime/explicit)
  operator T() const { return *static_cast<T*>(this); }  // upcast
};

}  // namespace math
}  // namespace stan
#endif
