#ifndef STAN_MATH_FWD_FUN_TCROSSPROD_HPP
#define STAN_MATH_FWD_FUN_TCROSSPROD_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/transpose.hpp>
#include <stan/math/fwd/fun/multiply.hpp>

namespace stan {
namespace math {

template <typename EigMat, require_eigen_vt<is_fvar, EigMat>* = nullptr>
inline auto tcrossprod(const EigMat& m) {
  using ret_type = plain_type_t<decltype(multiply(m, m.transpose()))>;
  if (m.rows() == 0) {
    return ret_type{};
  }
  const auto& m_ref = to_ref(m);
  return ret_type(multiply(m_ref, m_ref.transpose()));
}

}  // namespace math
}  // namespace stan
#endif
