#ifndef STAN_MATH_FWD_SCAL_FUN_TO_FVAR_HPP
#define STAN_MATH_FWD_SCAL_FUN_TO_FVAR_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {

template <typename T>
inline fvar<T> to_fvar(const T& x) {
  return fvar<T>(x);
}

template <typename T>
inline fvar<T> to_fvar(const fvar<T>& x) {
  return x;
}

template <typename T1, typename T2, enable_if_all_not_same<T2, T2>* = nullptr,
          enable_if_fvar<T2>* = nullptr>
inline fvar<T2> to_fvar(const T1& x, const T2& x2) {
  return fvar<T2>(x);
}

template <typename T1, typename T2, enable_if_all_same<T1, T2>* = nullptr>
inline T1& to_fvar(T1& x, const T2& x2) {
  return x;
}

template <typename T1, typename T2, enable_if_not_same<T1, T2>* = nullptr,
          enable_if_not_arithmetic<T2>* = nullptr>
inline fvar<T2> to_fvar(const fvar<T1>& x, const fvar<T2>& x2) {
  return fvar<T2>(x);
}

}  // namespace math
}  // namespace stan
#endif
