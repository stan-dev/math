#ifndef STAN_MATH_PRIM_FUN_POSITIVE_ORDERED_FREE_HPP
#define STAN_MATH_PRIM_FUN_POSITIVE_ORDERED_FREE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the vector of unconstrained scalars that transform to
 * the specified positive ordered vector.
 *
 * <p>This function inverts the constraining operation defined in
 * <code>positive_ordered_constrain(Matrix)</code>,
 *
 * @tparam T type of elements in the vector
 * @param y Vector of positive, ordered scalars.
 * @return Free vector that transforms into the input vector.
 * @throw std::domain_error if y is not a vector of positive,
 *   ordered scalars.
 */
template <typename EigVec, require_eigen_col_vector_t<EigVec>* = nullptr>
Eigen::Matrix<value_type_t<EigVec>, Eigen::Dynamic, 1> positive_ordered_free(
    const EigVec& y) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using std::log;
  using size_type = Eigen::Index;

  check_positive_ordered("stan::math::positive_ordered_free",
                         "Positive ordered variable", y);
  size_type k = y.size();
  Matrix<value_type_t<EigVec>, Dynamic, 1> x(k);
  if (k == 0) {
    return x;
  }
  x.coeffRef(0) = log(y.coeff(0));
  for (size_type i = 1; i < k; ++i) {
    x.coeffRef(i) = log(y.coeff(i) - y.coeff(i - 1));
  }
  return x;
}

}  // namespace math
}  // namespace stan

#endif
