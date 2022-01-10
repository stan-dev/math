#ifndef STAN_MATH_PRIM_PROB_EXPONENTIAL_QF_HPP
#define STAN_MATH_PRIM_PROB_EXPONENTIAL_QF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <vector>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * The quantile of an exponential density for p with the specified
 * inverse scale parameter.
 * Inverse scale parameter must be greater than 0.
 * p must be bounded by 0 and 1.
 *
 \f[
   q = -\log(1 - p) / \beta
 \f]
 *
 * @tparam Tp type of probability input
 * @tparam Tbeta type of inverse scale
 * @param p A scalar variable.
 * @param beta Inverse scale parameter.
 * @throw std::domain_error if beta is not greater than 0.
 * @throw std::domain_error if y is not greater than or equal to 0.
 */
template <typename Tp, typename Tbeta,
          require_all_arithmetic_t<Tp, Tbeta>* = nullptr>
inline auto exponential_qf(const Tp& p, const Tbeta& beta) {
  static const char* function = "exponential_qf";
  check_positive_finite(function, "Inverse scale parameter", beta);
  check_bounded(function, "Probability parameter", p, 0, 1);
  return -log1p(-p) / beta;
}

/** \ingroup prob_dists
 * The quantile of an exponential density for p with the specified
 * inverse scale parameter.
 * Inverse scale parameter must be greater than 0.
 * p must be bounded by 0 and 1.
 *
 * Specialisation for use where any input is an Eigen vector
 *
 * @tparam Tp type of probability input
 * @tparam Tbeta type of inverse scale
 * @param p A vector or scalar of probabilities.
 * @param beta Vector or scalar of inverse scale parameter.
 * @throw std::domain_error if beta is not greater than 0.
 * @throw std::domain_error if y is not greater than or equal to 0.
 */
template <typename Tp, typename Tbeta,
          require_all_st_arithmetic<Tp, Tbeta>* = nullptr,
          require_any_eigen_vector_t<Tp, Tbeta>* = nullptr>
inline auto exponential_qf(const Tp& p, const Tbeta& beta) {
  static const char* function = "exponential_qf";
  ref_type_t<Tp> p_ref = p;
  ref_type_t<Tbeta> beta_ref = beta;
  check_positive_finite(function, "Inverse scale parameter", beta_ref);
  check_bounded(function, "Probability parameter", p_ref, 0, 1);
  return (-log1p(-as_array_or_scalar(p_ref)) / as_array_or_scalar(beta_ref))
      .matrix();
}

}  // namespace math
}  // namespace stan

#endif
