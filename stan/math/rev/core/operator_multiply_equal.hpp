#ifndef STAN_MATH_REV_CORE_OPERATOR_MULTIPLY_EQUAL_HPP
#define STAN_MATH_REV_CORE_OPERATOR_MULTIPLY_EQUAL_HPP

#include <stan/math/rev/core/var.hpp>
#include <stan/math/rev/core/operator_multiplication.hpp>
#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {

template <typename T, typename VariType>
inline var_value<T, VariType>& var_value<T, VariType>::operator*=(const var_value<T, VariType>& b) {
  vi_ = new internal::multiply_vv_vari(vi_, b.vi_);
  return *this;
}

template <typename T, typename VariType>
inline var_value<T, VariType>& var_value<T, VariType>::operator*=(T b) {
  if (b == 1.0) {
    return *this;
  }
  vi_ = new internal::multiply_vd_vari(vi_, b);
  return *this;
}

}  // namespace math
}  // namespace stan
#endif
