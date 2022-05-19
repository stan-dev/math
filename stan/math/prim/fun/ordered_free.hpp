#ifndef STAN_MATH_PRIM_FUN_ORDERED_FREE_HPP
#define STAN_MATH_PRIM_FUN_ORDERED_FREE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the vector of unconstrained scalars that transform to
 * the specified positive ordered vector.
 *
 * <p>This function inverts the constraining operation defined in
 * <code>ordered_constrain(Matrix)</code>,
 *
 * @tparam T type of elements in the vector
 * @param y Vector of positive, ordered scalars.
 * @return Free vector that transforms into the input vector.
 * @throw std::domain_error if y is not a vector of positive,
 *   ordered scalars.
 */
template <typename EigVec, require_eigen_col_vector_t<EigVec>* = nullptr>
plain_type_t<EigVec> ordered_free(const EigVec& y) {
  const auto& y_ref = to_ref(y);
  check_ordered("stan::math::ordered_free", "Ordered variable", y_ref);
  using std::log;
  Eigen::Index k = y.size();
  plain_type_t<EigVec> x(k);
  if (k == 0) {
    return x;
  }
  x[0] = y_ref[0];
  for (Eigen::Index i = 1; i < k; ++i) {
    x.coeffRef(i) = log(y_ref.coeff(i) - y_ref.coeff(i - 1));
  }
  return x;
}

/**
 * Overload of `ordered_free()` to untransform each Eigen vector
 * in a standard vector.
 * @tparam T A standard vector with with a \ref stan::value_type which inherits
 * from `Eigen::MatrixBase` with compile time rows or columns equal to 1.
 * @param x The standard vector to untransform.
 */
template <typename T, require_std_vector_t<T>* = nullptr>
auto ordered_free(const T& x) {
  return apply_vector_unary<T>::apply(x,
                                      [](auto&& v) { return ordered_free(v); });
}

}  // namespace math
}  // namespace stan

#endif
