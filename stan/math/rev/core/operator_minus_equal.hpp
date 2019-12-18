#ifndef STAN_MATH_REV_CORE_OPERATOR_MINUS_EQUAL_HPP
#define STAN_MATH_REV_CORE_OPERATOR_MINUS_EQUAL_HPP

#include <stan/math/rev/core/var.hpp>
#include <stan/math/rev/core/operator_subtraction.hpp>
#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {

template <typename T, require_var_t<T>...>
inline var& var::operator-=(T&& b) {
  vi_ = new internal::subtract_vv_vari(vi_, b.vi_);
  return *this;
}

template <typename T, require_arithmetic_t<T>...>
inline var& var::operator-=(T b) {
  if (b == 0.0) {
    return *this;
  }
  vi_ = new internal::subtract_vd_vari(vi_, b);
  return *this;
}

}  // namespace math
}  // namespace stan
#endif
