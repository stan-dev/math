#ifndef STAN_MATH_REV_CORE_OPERATOR_MULTIPLY_EQUAL_HPP
#define STAN_MATH_REV_CORE_OPERATOR_MULTIPLY_EQUAL_HPP

#include <stan/math/rev/core/var.hpp>
#include <stan/math/rev/core/operator_multiplication.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>

namespace stan {
namespace math {

template <typename T>
inline var_value<T>& var_value<T, require_vt_floating_point<T>>::operator*=(
    const var_value<T>& b) {
  vi_ = new internal::multiply_vari<var_value<T>, var_value<T>>(vi_, b.vi_);
  return *this;
}

template <typename T>
template <typename Arith, require_vt_arithmetic<Arith>...>
inline var_value<T>& var_value<T, require_vt_floating_point<T>>::operator*=(
    const Arith& b) {
  vi_ = new internal::multiply_vari<var_value<T>, Arith>(vi_, b);
  return *this;
}

}  // namespace math
}  // namespace stan
#endif
