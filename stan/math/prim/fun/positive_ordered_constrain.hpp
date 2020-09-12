#ifndef STAN_MATH_PRIM_FUN_POSITIVE_ORDERED_CONSTRAIN_HPP
#define STAN_MATH_PRIM_FUN_POSITIVE_ORDERED_CONSTRAIN_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return an increasing positive ordered vector derived from the specified
 * free vector.  The returned constrained vector will have the
 * same dimensionality as the specified free vector.
 *
 * @tparam T type of elements in the vector
 * @param x Free vector of scalars.
 * @return Positive, increasing ordered vector.
 */
template <typename EigVec, require_eigen_col_vector_t<EigVec>* = nullptr>
Eigen::Matrix<value_type_t<EigVec>, Eigen::Dynamic, 1>
positive_ordered_constrain(const EigVec& x) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using std::exp;
  using size_type = Eigen::Index;

  size_type k = x.size();
  Matrix<value_type_t<EigVec>, Dynamic, 1> y(k);
  if (k == 0) {
    return y;
  }
  y.coeffRef(0) = exp(x.coeff(0));
  for (size_type i = 1; i < k; ++i) {
    y.coeffRef(i) = y.coeff(i - 1) + exp(x.coeff(i));
  }
  return y;
}

/**
 * Return a positive valued, increasing positive ordered vector derived
 * from the specified free vector and increment the specified log
 * probability reference with the log absolute Jacobian determinant
 * of the transform.  The returned constrained vector
 * will have the same dimensionality as the specified free vector.
 *
 * @tparam T type of elements in the vector
 * @param x Free vector of scalars.
 * @param lp Log probability reference.
 * @return Positive, increasing ordered vector.
 */
template <typename EigVec, require_eigen_col_vector_t<EigVec>* = nullptr>
inline Eigen::Matrix<value_type_t<EigVec>, Eigen::Dynamic, 1>
positive_ordered_constrain(const EigVec& x, value_type_t<EigVec>& lp) {
  lp += sum(x);
  return positive_ordered_constrain(x);
}

}  // namespace math
}  // namespace stan

#endif
