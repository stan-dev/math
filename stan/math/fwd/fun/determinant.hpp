#ifndef STAN_MATH_FWD_FUN_DETERMINANT_HPP
#define STAN_MATH_FWD_FUN_DETERMINANT_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/fwd/core.hpp>

namespace stan {
namespace math {

template <typename T, require_eigen_vt<is_fvar, T>* = nullptr>
inline value_type_t<T> determinant(const T& m) {
  check_square("determinant", "m", m);

  const typename value_type_t<T>::Scalar vals = m.val().determinant();
  return value_type_t<T>(vals, vals * (m.val().inverse() * m.d()).trace());
}

}  // namespace math
}  // namespace stan
#endif
