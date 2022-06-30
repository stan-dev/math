#ifndef STAN_MATH_PRIM_PROB_POISSON_GAMMA_LCDF_HPP
#define STAN_MATH_PRIM_PROB_POISSON_GAMMA_LCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/elt_divide.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/prob/neg_binomial_2_lcdf.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Returns the log of the cumulative distribution function (LCDF) for a Poisson
 * random variate with a Gamma prior for the rate parameter
 *
 * @tparam T_n Type of count outcome
 * @tparam T_shape Type of shape parameter for the Gamma prior
 * @tparam T_inv_scale Type of rate parameter for the Gamma prior
 * @param n Discrete count outcome
 * @param alpha Shape parameter for the Gamma prior
 * @param beta Rate parameter for the Gamma prior
 * @return Log CDF
 */
template <typename T_n, typename T_shape, typename T_inv_scale>
return_type_t<T_shape, T_inv_scale> poisson_gamma_lcdf(const T_n& n,
                                                      const T_shape& alpha,
                                                      const T_inv_scale& beta) {
  static const char* function = "poisson_gamma_lcdf";
  const auto& n_ref = to_ref(n);
  const auto& alpha_ref = to_ref(as_column_vector_or_scalar(alpha));
  const auto& beta_ref = to_ref(as_column_vector_or_scalar(beta));

  check_nonnegative(function, "Random variable", n_ref);
  check_positive_finite(function, "Shape parameter", alpha_ref);
  check_positive_finite(function, "Inverse scale parameter", beta_ref);

  return neg_binomial_2_lcdf(n_ref, elt_divide(alpha_ref, beta_ref),
                             alpha_ref);
}

}  // namespace math
}  // namespace stan
#endif
