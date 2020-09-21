#ifndef STAN_MATH_REV_CORE_OPERATOR_MULTIPLY_EQUAL_HPP
#define STAN_MATH_REV_CORE_OPERATOR_MULTIPLY_EQUAL_HPP

#include <stan/math/rev/core/var.hpp>
#include <stan/math/rev/core/operator_multiplication.hpp>
#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {

template <typename T>
inline var_value<T>& var_value<T, require_floating_point_t<T>>::operator*=(
    const var_value<T>& b) {
  vi_ = new internal::multiply_vv_vari(vi_, b.vi_);
  return *this;
}

template <typename T>
inline var_value<T>& var_value<T, require_floating_point_t<T>>::operator*=(
    T b) {
  if (b == 1.0) {
    return *this;
  }
  vi_ = new internal::multiply_vd_vari(vi_, b);
  return *this;
}

}  // namespace math
}  // namespace stan
#endif
