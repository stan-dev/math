#ifndef STAN_MATH_PRIM_FUN_UNIT_VECTOR_CONSTRAIN_HPP
#define STAN_MATH_PRIM_FUN_UNIT_VECTOR_CONSTRAIN_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/dot_self.hpp>
#include <stan/math/prim/fun/sqrt.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the unit length vector corresponding to the free vector y.
 *
 * See <a
 * href="https://en.wikipedia.org/wiki/N-sphere#Generating_random_points">the
 * Wikipedia page on generating random points on an N-sphere</a>.
 *
 * @tparam EigMat type inheriting from `EigenBase` that does not have an fvar
 *  scalar type.
 * @param y vector of K unrestricted variables
 * @return Unit length vector of dimension K
 */
template <typename EigMat, require_eigen_t<EigMat>* = nullptr,
          require_not_vt_fvar<EigMat>* = nullptr>
inline auto unit_vector_constrain(const EigMat& y) {
  using std::sqrt;
  const auto& y_ref = to_ref(y);
  check_vector("unit_vector_constrain", "y", y_ref);
  check_nonzero_size("unit_vector_constrain", "y", y_ref);
  value_type_t<EigMat> SN = dot_self(y_ref);
  check_positive_finite("unit_vector_constrain", "norm", SN);
  return (y_ref / sqrt(SN)).eval();
}

/**
 * Return the unit length vector corresponding to the free vector y.
 * See https://en.wikipedia.org/wiki/N-sphere#Generating_random_points
 *
 * @tparam EigMat type inheriting from `EigenBase` that does not have an fvar
 *  scalar type.
 *
 * @param y vector of K unrestricted variables
 * @return Unit length vector of dimension K
 * @param lp Log probability reference to increment.
 */
template <typename EigMat, typename LP, require_eigen_t<EigMat>* = nullptr,
          require_not_vt_fvar<EigMat>* = nullptr>
inline auto unit_vector_constrain(const EigMat& y, LP& lp) {
  const auto& y_ref = to_ref(y);
  using std::sqrt;
  check_vector("unit_vector_constrain", "y", y_ref);
  check_nonzero_size("unit_vector_constrain", "y", y_ref);
  value_type_t<EigMat> SN = dot_self(y_ref);
  check_positive_finite("unit_vector_constrain", "norm", SN);
  lp -= 0.5 * SN;
  return (y_ref / sqrt(SN)).eval();
}

}  // namespace math
}  // namespace stan

#endif
