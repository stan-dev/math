#ifndef STAN_MATH_REV_CORE_OPERATOR_DIVIDE_EQUAL_HPP
#define STAN_MATH_REV_CORE_OPERATOR_DIVIDE_EQUAL_HPP

#include <stan/math/rev/core/var.hpp>
#include <stan/math/rev/core/operator_division.hpp>
#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {
template <typename T, typename VariType>
template <typename OtherT, typename OtherVariType>
inline var_value<T, VariType>& var_value<T, VariType>::operator/=(
    const var_value<OtherT, OtherVariType>& b) {
  vi_ = new internal::divide_vv_vari(vi_, b.vi_);
  return *this;
}

template <typename T, typename VariType>
template <typename S, require_convertible_t<S&, T>*>
inline var_value<T, VariType>& var_value<T, VariType>::operator/=(
    const S& b) {
  if (b == 1.0) {
    return *this;
  }
  vi_ = new internal::divide_vd_vari(vi_, b);
  return *this;
}

}  // namespace math
}  // namespace stan
#endif
