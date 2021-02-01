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
template <typename T, require_eigen_col_vector_t<T>* = nullptr,
          require_not_vt_autodiff<T>* = nullptr>
inline auto unit_vector_constrain(const T& y) {
  using std::sqrt;
  check_nonzero_size("unit_vector_constrain", "y", y);
  return make_holder(
      [](const auto& y_ref) {
        value_type_t<T> SN = dot_self(y_ref);
        check_positive_finite("unit_vector_constrain", "norm", SN);
        return y_ref / sqrt(SN);
      },
      to_ref(y));
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
template <typename T1, typename T2, require_eigen_col_vector_t<T1>* = nullptr,
          require_all_not_vt_autodiff<T1, T2>* = nullptr>
inline plain_type_t<T1> unit_vector_constrain(const T1& y, T2& lp) {
  using std::sqrt;
  check_nonzero_size("unit_vector_constrain", "y", y);
  return make_holder(
      [](const auto& y_ref, auto& lp) {
        value_type_t<T1> SN = dot_self(y_ref);
        check_positive_finite("unit_vector_constrain", "norm", SN);
        lp -= 0.5 * SN;
        return y_ref / sqrt(SN);
      },
      to_ref(y), lp);
}

}  // namespace math
}  // namespace stan

#endif
