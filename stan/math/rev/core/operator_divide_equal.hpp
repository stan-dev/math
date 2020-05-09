#ifndef STAN_MATH_REV_CORE_OPERATOR_DIVIDE_EQUAL_HPP
#define STAN_MATH_REV_CORE_OPERATOR_DIVIDE_EQUAL_HPP

#include <stan/math/rev/core/var.hpp>
#include <stan/math/rev/core/operator_division.hpp>
#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {

template <typename T>
inline var_value<T>& var_value<T>::operator/=(const var_value<T>& b) {
  vi_ = new internal::divide_vari<T, vari_value<T>, vari_value<T>>(vi_, b.vi_);
  return *this;
}

template <typename T>
template <typename Arith, require_arithmetic_t<Arith>...>
inline var_value<T>& var_value<T>::operator/=(const Arith& b) {
  if (b == 1.0) {
    return *this;
  }
  vi_ = new internal::divide_vari<T, vari_value<T>, Arith>(vi_, b);
  return *this;
}

}  // namespace math
}  // namespace stan
#endif
