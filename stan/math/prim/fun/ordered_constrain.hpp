#ifndef STAN_MATH_PRIM_FUN_ORDERED_CONSTRAIN_HPP
#define STAN_MATH_PRIM_FUN_ORDERED_CONSTRAIN_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return an increasing ordered vector derived from the specified
 * free vector.  The returned constrained vector will have the
 * same dimensionality as the specified free vector.
 *
 * @tparam T type of the vector
 * @param x Free vector of scalars.
 * @return Positive, increasing ordered vector.
 * @tparam T Type of scalar.
 */
template <typename EigVec, require_eigen_col_vector_t<EigVec>* = nullptr,
          require_not_st_var<EigVec>* = nullptr>
plain_type_t<EigVec> ordered_constrain(const EigVec& x) {
  using std::exp;
  Eigen::Index k = x.size();
  plain_type_t<EigVec> y(k);
  const auto& x_ref = to_ref(x);
  if (unlikely(k == 0)) {
    return y;
  }
  y[0] = x_ref[0];
  for (Eigen::Index i = 1; i < k; ++i) {
    y.coeffRef(i) = y.coeff(i - 1) + exp(x_ref.coeff(i));
  }
  return y;
}

/**
 * Return a positive valued, increasing ordered vector derived
 * from the specified free vector and increment the specified log
 * probability reference with the log absolute Jacobian determinant
 * of the transform.  The returned constrained vector
 * will have the same dimensionality as the specified free vector.
 *
 * @tparam T type of the vector
 * @param x Free vector of scalars.
 * @param lp Log probability reference.
 * @return Positive, increasing ordered vector.
 */
template <typename EigVec, require_eigen_col_vector_t<EigVec>* = nullptr>
auto ordered_constrain(const EigVec& x, value_type_t<EigVec>& lp) {
  const auto& x_ref = to_ref(x);
  if (likely(x.size() > 1)) {
    lp += x_ref.tail(x.size() - 1).sum();
  }
  return ordered_constrain(x_ref);
}

}  // namespace math
}  // namespace stan

#endif
