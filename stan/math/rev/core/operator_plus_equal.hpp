#ifndef STAN_MATH_REV_CORE_OPERATOR_PLUS_EQUAL_HPP
#define STAN_MATH_REV_CORE_OPERATOR_PLUS_EQUAL_HPP

#include <stan/math/rev/core/var.hpp>
#include <stan/math/rev/core/operator_addition.hpp>
#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {

template <typename T>
inline var_type<T>& var_type<T>::operator+=(var_type<T> b) {
  vi_ = new internal::add_vari<T, vari_type<T>, vari_type<T>>(vi_, b.vi_);
  return *this;
}

template <typename T>
template <typename Arith, require_vt_arithmetic<Arith>...>
inline var_type<T>& var_type<T>::operator+=(Arith b) {
  // TODO: No internal from elsewhere!
  if (internal::is_any_equal(b, 0.0)) {
    return *this;
  }
  vi_ = new internal::add_vari<double, vari_type<T>, Arith>(vi_, b);
  return *this;
}

}  // namespace math
}  // namespace stan
#endif
