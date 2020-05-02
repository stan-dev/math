#ifndef STAN_MATH_FWD_FUN_DETERMINANT_HPP
#define STAN_MATH_FWD_FUN_DETERMINANT_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/fwd/core.hpp>

namespace stan {
namespace math {

template <typename EigMat, require_eigen_vt<is_fvar, EigMat>* = nullptr>
inline value_type_t<EigMat> determinant(const EigMat& m) {
  check_square("determinant", "m", m);

  const typename value_type_t<EigMat>::Scalar vals = m.val().determinant();
  return {vals, vals * (m.val().inverse() * m.d()).trace()};
}

}  // namespace math
}  // namespace stan
#endif
