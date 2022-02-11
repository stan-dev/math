#ifndef STAN_MATH_OPENCL_PRIM_OFFSET_MULTIPLIER_CONSTRAIN_HPP
#define STAN_MATH_OPENCL_PRIM_OFFSET_MULTIPLIER_CONSTRAIN_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/prim/sum.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/kernel_generator.hpp>

namespace stan {
namespace math {

/**
 * Return the linearly transformed value for the specified unconstrained input
 * and specified offset and multiplier.
 *
 * <p>The transform applied is
 *
 * <p>\f$f(x) = mu + sigma * x\f$
 *
 * <p>where mu is the offset and sigma is the multiplier.
 *
 * <p>If the offset is zero and the multiplier is one this
 * reduces to <code>identity_constrain(x)</code>.
 *
 * @tparam T type of unconstrained input
 * @tparam M type of offset
 * @tparam S type of multiplier
 * @param[in] x Unconstrained input
 * @param[in] mu offset of constrained output
 * @param[in] sigma multiplier of constrained output
 * @return linear transformed value corresponding to inputs
 * @throw std::domain_error if sigma <= 0
 * @throw std::domain_error if mu is not finite
 */
template <typename T, typename M, typename S,
          require_all_kernel_expressions_t<T, M, S>* = nullptr,
          require_any_not_stan_scalar_t<T, M, S>* = nullptr>
inline auto offset_multiplier_constrain(const T& x, const M& mu,
                                        const S& sigma) {
  using std::isfinite;
  const char* function = "offset_multiplier_constrain(OpenCL)";
  check_consistent_sizes(function, "offset", mu, "multiplier", sigma,
                         "parameter", x);
  auto check_mu = check_cl(function, "offset", mu, "finite");
  auto check_sigma = check_cl(function, "multiplier", sigma, "positive finite");
  matrix_cl<double> res;
  results(check_mu, check_sigma, res) = expressions(
      isfinite(mu), isfinite(sigma) && sigma > 0, elt_multiply(x, sigma) + mu);
  return res;
}

/**
 * Return the linearly transformed value for the specified unconstrained input
 * and specified offset and multiplier.
 *
 * <p>The transform applied is
 *
 * <p>\f$f(x) = mu + sigma * x\f$
 *
 * <p>where mu is the offset and sigma is the multiplier.
 *
 * <p>If the offset is zero and the multiplier is one this
 * reduces to <code>identity_constrain(x)</code>.
 *
 * @tparam T type of unconstrained input
 * @tparam M type of offset
 * @tparam S type of multiplier
 * @param[in] x Unconstrained input
 * @param[in] mu offset of constrained output
 * @param[in] sigma multiplier of constrained output
 * @param[in,out] lp Reference to log probability to increment.
 * @return linear transformed value corresponding to inputs
 * @throw std::domain_error if sigma <= 0
 * @throw std::domain_error if mu is not finite
 */
template <typename T, typename M, typename S,
          require_all_kernel_expressions_t<T, M, S>* = nullptr,
          require_any_not_stan_scalar_t<T, M, S>* = nullptr>
inline auto offset_multiplier_constrain(const T& x, const M& mu, const S& sigma,
                                        double& lp) {
  using std::isfinite;
  const char* function = "offset_multiplier_constrain(OpenCL)";
  check_consistent_sizes(function, "offset", mu, "multiplier", sigma,
                         "parameter", x);
  auto check_mu = check_cl(function, "offset", mu, "finite");
  auto check_sigma = check_cl(function, "multiplier", sigma, "positive finite");
  matrix_cl<double> res;
  matrix_cl<double> lp_inc;
  results(check_mu, check_sigma, res, lp_inc)
      = expressions(isfinite(mu), isfinite(sigma) && sigma > 0,
                    elt_multiply(x, sigma) + mu, sum_2d(log(sigma)));
  lp += sum(from_matrix_cl(lp_inc));
  return res;
}

}  // namespace math
}  // namespace stan
#endif
#endif
