#ifndef STAN_MATH_PRIM_PROB_INV_GAUSSIAN_LCCDF_HPP
#define STAN_MATH_PRIM_PROB_INV_GAUSSIAN_LCCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/log_diff_exp.hpp>
#include <stan/math/prim/fun/select.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/sqrt.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>
#include <stan/math/prim/functor/apply_scalar_binary.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <stan/math/prim/prob/inv_gaussian_lcdf.hpp>
#include <cmath>
#include <utility>

namespace stan {
namespace math {
namespace internal {

/**
 * Return log_diff_exp(a, b), or negative infinity when b >= a. The survivor
 * is a difference of two terms that meet in the far upper tail, where
 * rounding can invert their order and zero is the correct answer.
 */
template <typename T1, typename T2,
          require_all_stan_scalar_t<T1, T2>* = nullptr>
inline return_type_t<T1, T2> log_diff_exp_guarded(const T1& a, const T2& b) {
  if (unlikely(value_of_rec(b) >= value_of_rec(a))) {
    return NEGATIVE_INFTY;
  }
  return log_diff_exp(a, b);
}

/**
 * A vectorized version of log_diff_exp_guarded().
 */
template <typename T1, typename T2, require_any_container_t<T1, T2>* = nullptr>
inline auto log_diff_exp_guarded(T1&& a, T2&& b) {
  return apply_scalar_binary(
      [](auto&& c, auto&& d) {
        return log_diff_exp_guarded(std::forward<decltype(c)>(c),
                                    std::forward<decltype(d)>(d));
      },
      std::forward<T1>(a), std::forward<T2>(b));
}

}  // namespace internal

/** \ingroup prob_dists
 * Returns the inverse Gaussian log complementary cumulative distribution
 * function for the given random variable, mean and shape. Given containers
 * of matching sizes, returns the log of the product of complementary
 * probabilities.
 *
 * <p>The survivor function is \f$S(y \mid \mu, \lambda) = \Phi(-z_1)
 * - e^{2\lambda/\mu}\Phi(-z_2)\f$ with
 * \f$z_{1,2} = \sqrt{\lambda/y}\,(y/\mu \mp 1)\f$. Both terms are evaluated
 * as lower-tail log-CDFs, which stay accurate far into the tail, and the
 * difference is taken in log space. See <code>inv_gaussian_lcdf</code>.
 *
 * <p>Deep in the upper tail the two terms meet within the spacing of a
 * double and the guarded difference returns \f$-\infty\f$. Over
 * \f$\mu \in [10^{-3}, 10^3]\f$, \f$\lambda \in [10^{-3}, 10^{16}]\f$
 * that first happens at \f$\log S = -5 \times 10^5\f$, so the survivor is
 * zero to any representable precision there.
 *
 * <p>Both boundaries are handled elementwise: \f$y = 0\f$ contributes
 * \f$\log 1 = 0\f$ and \f$y = \infty\f$ contributes \f$-\infty\f$, so a
 * container mixing a boundary with ordinary observations sums to the same
 * value as its elements taken one at a time. The partials are zero at both,
 * and wherever the survivor has underflowed.
 *
 * @tparam T_y type of scalar
 * @tparam T_loc type of mean parameter
 * @tparam T_shape type of shape parameter
 * @param y scalar or container of scalars
 * @param mu mean parameter
 * @param lambda shape parameter
 * @return the log of the product of complementary probabilities
 * @throw std::domain_error if the mean or the shape is not positive and
 * finite, if the random variable is negative, or if any argument is NaN.
 * @throw std::invalid_argument if container sizes mismatch.
 */
template <typename T_y, typename T_loc, typename T_shape,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_loc, T_shape>* = nullptr>
inline return_type_t<T_y, T_loc, T_shape> inv_gaussian_lccdf(
    const T_y& y, const T_loc& mu, const T_shape& lambda) {
  using T_partials_return = partials_return_t<T_y, T_loc, T_shape>;
  using T_y_ref = ref_type_if_not_constant_t<T_y>;
  using T_mu_ref = ref_type_if_not_constant_t<T_loc>;
  using T_lambda_ref = ref_type_if_not_constant_t<T_shape>;
  static constexpr const char* function = "inv_gaussian_lccdf";
  check_consistent_sizes(function, "Random variable", y, "Mean parameter", mu,
                         "Shape parameter", lambda);

  T_y_ref y_ref = y;
  T_mu_ref mu_ref = mu;
  T_lambda_ref lambda_ref = lambda;

  decltype(auto) y_val = to_ref(as_value_column_array_or_scalar(y_ref));
  decltype(auto) mu_val = to_ref(as_value_column_array_or_scalar(mu_ref));
  decltype(auto) lambda_val
      = to_ref(as_value_column_array_or_scalar(lambda_ref));

  check_nonnegative(function, "Random variable", y_val);
  check_positive_finite(function, "Mean parameter", mu_val);
  check_positive_finite(function, "Shape parameter", lambda_val);

  if (size_zero(y, mu, lambda)) {
    return 0;
  }

  auto ops_partials = make_partials_propagator(y_ref, mu_ref, lambda_ref);

  // Boundary masks; see inv_gaussian_lcdf.
  const auto& is_inf = to_ref(y_val == INFTY);
  const auto& is_bdry = to_ref((y_val == 0) || (y_val == INFTY));

  constexpr bool any_ad = is_any_autodiff_v<T_y, T_loc, T_shape>;

  const auto& inv_mu = to_ref_if<any_ad>(inv(mu_val));
  const auto& inv_y = to_ref_if<any_ad>(inv(y_val));
  const auto& sqrt_lambda_over_y = to_ref_if<any_ad>(sqrt(lambda_val * inv_y));
  const auto& y_over_mu = to_ref_if<any_ad>(y_val * inv_mu);
  const auto& z1 = to_ref(sqrt_lambda_over_y * (y_over_mu - 1.0));
  const auto& z2 = to_ref(sqrt_lambda_over_y * (y_over_mu + 1.0));

  const auto& log_upper = to_ref(internal::log_scaled_upper_term(z1, z2));
  const auto& lccdf_elt = to_ref(select(
      is_inf, T_partials_return(NEGATIVE_INFTY),
      internal::log_diff_exp_guarded(internal::log_Phi(-z1), log_upper)));

  T_partials_return lccdf = sum(lccdf_elt);

  if constexpr (any_ad) {
    // Same phi collapse as in the lcdf; the partials are its sign flip, with
    // the survivor in place of the CDF.
    //
    // Each partial is masked separately: at y == inf the shape partial is
    // 0 / 0 on its own. The inner select covers elements whose survivor has
    // underflowed to -inf.
    const auto& is_underflow = to_ref(lccdf_elt == NEGATIVE_INFTY);
    const auto& w_dens = to_ref(exp(internal::log_phi(z1) - lccdf_elt));
    const auto& w_upper = to_ref(exp(log_upper - lccdf_elt));
    if constexpr (is_autodiff_v<T_y>) {
      partials<0>(ops_partials)
          = select(is_bdry, T_partials_return(0),
                   select(is_underflow, T_partials_return(0),
                          -w_dens * sqrt_lambda_over_y * inv_y));
    }
    if constexpr (is_autodiff_v<T_loc>) {
      partials<1>(ops_partials)
          = select(is_bdry, T_partials_return(0),
                   select(is_underflow, T_partials_return(0),
                          2.0 * lambda_val * w_upper * square(inv_mu)));
    }
    if constexpr (is_autodiff_v<T_shape>) {
      partials<2>(ops_partials) = select(
          is_bdry, T_partials_return(0),
          select(is_underflow, T_partials_return(0),
                 w_dens * inv_y / sqrt_lambda_over_y - 2.0 * w_upper * inv_mu));
    }
  }
  return ops_partials.build(lccdf);
}

}  // namespace math
}  // namespace stan
#endif
