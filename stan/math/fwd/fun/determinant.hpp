#ifndef STAN_MATH_FWD_FUN_DETERMINANT_HPP
#define STAN_MATH_FWD_FUN_DETERMINANT_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/to_ref.hpp>

namespace stan {
namespace math {

template <typename EigMat, require_eigen_vt<is_fvar, EigMat>* = nullptr>
inline value_type_t<EigMat> determinant(const EigMat& m) {
  check_square("determinant", "m", m);
  const auto& m_ref = to_ref(m);

  const typename value_type_t<EigMat>::Scalar vals = m_ref.val().determinant();
  return {vals, vals * (m_ref.val().inverse() * m_ref.d()).trace()};
}

}  // namespace math
}  // namespace stan
#endif
