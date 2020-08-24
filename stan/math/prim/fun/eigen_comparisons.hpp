#ifndef STAN_MATH_PRIM_FUN_EIGEN_COMPARISONS_HPP
#define STAN_MATH_PRIM_FUN_EIGEN_COMPARISONS_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/value_of.hpp>

namespace stan {
namespace math {

template <typename T_a, typename T_b,
          require_any_eigen_vt<is_autodiff, T_a, T_b>* = nullptr,
          require_not_vt_same<T_a, T_b>* = nullptr>
auto operator<(const T_a& a, const T_b& b) {
  return operator<(value_of(a), value_of(b));
}

}  // namespace math
}  // namespace stan

#endif
