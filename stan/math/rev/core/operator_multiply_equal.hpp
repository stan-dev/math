#ifndef STAN_MATH_REV_CORE_OPERATOR_MULTIPLY_EQUAL_HPP
#define STAN_MATH_REV_CORE_OPERATOR_MULTIPLY_EQUAL_HPP

#include <stan/math/rev/core/var.hpp>
#include <stan/math/rev/core/operator_multiplication.hpp>
#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {

template <typename T>
inline var_type<T>& var_type<T>::operator*=(const var_type<T>& b) {
  vi_ = new internal::multiply_vari<T, vari_type<T>, vari_type<T>>(this->vi_, b.vi_);
  return *this;
}

template <typename T>
template <typename Arith, require_vt_arithmetic<Arith>...>
inline var_type<T>& var_type<T>::operator*=(Arith b) {
  if (internal::is_any_equal(b, 1.0)) {
    return *this;
  }
  vi_ = new internal::multiply_vari<T, vari_type<T>, Arith>(vi_, b);
  return *this;
}

}  // namespace math
}  // namespace stan
#endif
