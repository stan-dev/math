#ifndef STAN_MATH_PRIM_PROB_WEIBULL_LCDF_HPP
#define STAN_MATH_PRIM_PROB_WEIBULL_LCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/operands_and_partials.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Returns the Weibull log cumulative distribution function for the given
 * location and scale. Given containers of matching sizes, returns the
 * log sum of probabilities.
 *
 * @tparam T_y type of real parameter
 * @tparam T_shape type of shape parameter
 * @tparam T_scale type of scale parameter
 * @param y real parameter
 * @param alpha shape parameter
 * @param sigma scale parameter
 * @return log probability or log sum of probabilities
 * @throw std::domain_error if y is negative, alpha sigma is nonpositive
 */
template <typename T_y, typename T_shape, typename T_scale,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_shape, T_scale>* = nullptr>
return_type_t<T_y, T_shape, T_scale> weibull_lcdf(const T_y& y,
                                                  const T_shape& alpha,
                                                  const T_scale& sigma) {
  using T_partials_return = partials_return_t<T_y, T_shape, T_scale>;
  using T_y_ref = ref_type_t<T_y>;
  using T_alpha_ref = ref_type_t<T_shape>;
  using T_sigma_ref = ref_type_t<T_scale>;
  using std::pow;
  static const char* function = "weibull_lcdf";

  T_y_ref y_ref = y;
  T_alpha_ref alpha_ref = alpha;
  T_sigma_ref sigma_ref = sigma;

  check_nonnegative(function, "Random variable", y_ref);
  check_positive_finite(function, "Shape parameter", alpha_ref);
  check_positive_finite(function, "Scale parameter", sigma_ref);

  if (size_zero(y, alpha, sigma)) {
    return 0.0;
  }

  T_partials_return cdf_log(0);
  operands_and_partials<T_y_ref, T_alpha_ref, T_sigma_ref> ops_partials(
      y_ref, alpha_ref, sigma_ref);

  scalar_seq_view<T_y_ref> y_vec(y_ref);
  scalar_seq_view<T_alpha_ref> alpha_vec(alpha_ref);
  scalar_seq_view<T_sigma_ref> sigma_vec(sigma_ref);

  size_t max_size_seq_view = max_size(y, alpha, sigma);
  size_t size_y = stan::math::size(y);

  for (size_t i = 0; i < size_y; i++) {
    if (y_vec[i] == 0) {
      return LOG_ZERO;
    }
  }

  for (size_t i = 0; i < max_size_seq_view; i++) {
    const T_partials_return pow_n
        = pow(y_vec.val(i) / sigma_vec.val(i), alpha_vec.val(i));
    const T_partials_return exp_n = exp(-pow_n);
    const T_partials_return log1m_exp_n = log1m_exp(-pow_n);
    const T_partials_return rep_deriv = exp_n * pow_n / exp(log1m_exp_n);

    cdf_log += log1m_exp_n;

    if (!is_constant_all<T_y, T_scale, T_shape>::value) {
      const T_partials_return deriv_y_sigma = rep_deriv * alpha_vec.val(i);

      if (!is_constant_all<T_y>::value) {
        ops_partials.edge1_.partials_[i] += deriv_y_sigma / y_vec.val(i);
      }
      if (!is_constant_all<T_scale>::value) {
        ops_partials.edge3_.partials_[i] += -deriv_y_sigma / sigma_vec.val(i);
      }
    }
    if (!is_constant_all<T_shape>::value) {
      ops_partials.edge2_.partials_[i]
          += rep_deriv * log(y_vec.val(i) / sigma_vec.val(i));
    }
  }
  return ops_partials.build(cdf_log);
}
}  // namespace math
}  // namespace stan
#endif
