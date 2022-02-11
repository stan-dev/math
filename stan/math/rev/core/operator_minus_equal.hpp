#ifndef STAN_MATH_REV_CORE_OPERATOR_MINUS_EQUAL_HPP
#define STAN_MATH_REV_CORE_OPERATOR_MINUS_EQUAL_HPP

#include <stan/math/rev/core/var.hpp>
#include <stan/math/rev/core/operator_subtraction.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/size.hpp>
namespace stan {
namespace math {

template <typename T>
inline var_value<T>& var_value<T, require_floating_point_t<T>>::operator-=(
    const var_value<T>& b) {
  vi_ = (*this - b).vi_;
  return *this;
}

template <typename T>
inline var_value<T>& var_value<T, require_floating_point_t<T>>::operator-=(
    T b) {
  if (unlikely(b == 0.0)) {
    return *this;
  }
  vi_ = (*this - b).vi_;
  return *this;
}

}  // namespace math
}  // namespace stan
#endif
