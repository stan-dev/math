#ifndef STAN_MATH_PRIM_FUN_SIMPLEX_CONSTRAIN_HPP
#define STAN_MATH_PRIM_FUN_SIMPLEX_CONSTRAIN_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/inv_logit.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log1p_exp.hpp>
#include <stan/math/prim/fun/logit.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the simplex corresponding to the specified free vector.
 * A simplex is a vector containing values greater than or equal
 * to 0 that sum to 1.  A vector with (K-1) unconstrained values
 * will produce a simplex of size K.
 *
 * The transform is based on a centered stick-breaking process.
 *
 * @tparam Vec type deriving from `Eigen::MatgrixBase` with rows or columns
 * equal to 1.
 * @param y Free vector input of dimensionality K - 1.
 * @return Simplex of dimensionality K.
 */
template <typename Vec, require_eigen_vector_t<Vec>* = nullptr>
auto simplex_constrain(Vec&& y) {
  // cut & paste simplex_constrain(Eigen::Matrix, T) w/o Jacobian
  using std::log;
  const auto Km1 = y.size();
  const Eigen::Ref<const plain_type_t<Vec>>& y_ref = y;
  plain_type_t<Vec> x(Km1 + 1);
  value_type_t<Vec> stick_len(1.0);
  for (int k = 0; k < Km1; ++k) {
    const auto z_k = inv_logit(y_ref[k] - log(Km1 - k));
    x[k] = stick_len * z_k;
    stick_len -= x[k];
  }
  x[Km1] = stick_len;
  return x;
}

/**
 * Return the simplex corresponding to the specified free vector
 * and increment the specified log probability reference with
 * the log absolute Jacobian determinant of the transform.
 *
 * The simplex transform is defined through a centered
 * stick-breaking process.
 *
 * @tparam Vec type deriving from `Eigen::MatrixBase` with rows or columns
 * equal to 1.
 * @tparam T type of log probability.
 * @param y Free vector input of dimensionality K - 1.
 * @param lp Log probability reference to increment.
 * @return Simplex of dimensionality K.
 */
template <typename Vec, typename T, require_eigen_vector_t<Vec>* = nullptr>
auto simplex_constrain(Vec&& y, T& lp) {
  using std::log;
  const auto Km1 = y.size();  // K = Km1 + 1
  const Eigen::Ref<const plain_type_t<Vec>>& y_ref = y;
  plain_type_t<Vec> x(Km1 + 1);
  value_type_t<Vec> stick_len(1.0);
  for (int k = 0; k < Km1; ++k) {
    const auto eq_share = -log(Km1 - k);  // = logit(1.0/(Km1 + 1 - k));
    const auto adj_y_k = y_ref[k] + eq_share;
    const auto z_k = inv_logit(adj_y_k);
    x[k] = stick_len * z_k;
    lp += log(stick_len);
    lp -= log1p_exp(-adj_y_k);
    lp -= log1p_exp(adj_y_k);
    stick_len -= x[k];  // equivalently *= (1 - z_k);
  }
  x[Km1] = stick_len;  // no Jacobian contrib for last dim
  return x;
}

}  // namespace math
}  // namespace stan

#endif
