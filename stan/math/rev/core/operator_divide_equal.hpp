#ifndef STAN_MATH_REV_CORE_OPERATOR_DIVIDE_EQUAL_HPP
#define STAN_MATH_REV_CORE_OPERATOR_DIVIDE_EQUAL_HPP

#include <stan/math/rev/core/var.hpp>
#include <stan/math/rev/core/operator_division.hpp>
#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {
template <typename T>
inline var_type<T>& var_type<T>::operator/=(const var_type<T>& b) {
  vi_ = new internal::divide_vari<T, vari_type<T>, vari_type<T>>(vi_, b.vi_);
  return *this;
}

template <>
template <typename Arith, require_arithmetic_t<Arith>...>
inline var& var::operator/=(Arith b) {
  if (b == 1.0) {
    return *this;
  }
  vi_ = new internal::divide_vd_vari(vi_, b);
  return *this;
}

}  // namespace math
}  // namespace stan
#endif
