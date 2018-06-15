#ifndef STAN_MATH_PRIM_SCAL_META_RM_COMPLEX_HPP
#define STAN_MATH_PRIM_SCAL_META_RM_COMPLEX_HPP

#include <stan/math/prim/scal/meta/rm_zeroing.hpp>
#include <complex>

namespace stan {
namespace math {

template <class T>
struct complex;  // forward declare stan's complex

}  // namespace math

/// trait to remove the complex wrapper around a type
template <class T>
struct rm_complex {
  typedef rm_zeroing_t<T> type;
};
template <class T>
struct rm_complex<std::complex<T>> {
  typedef rm_zeroing_t<T> type;
};
template <class T>
struct rm_complex<stan::math::complex<T>> {
  typedef rm_zeroing_t<T> type;
};
template <class T>
using rm_complex_t = typename rm_complex<T>::type;

}  // namespace stan
#endif
