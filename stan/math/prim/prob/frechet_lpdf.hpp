#ifndef STAN_MATH_PRIM_PROB_FRECHET_LPDF_HPP
#define STAN_MATH_PRIM_PROB_FRECHET_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log1m.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/multiply_log.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <boost/random/weibull_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <stan/math/prim/functor/operands_and_partials.hpp>
#include <cmath>

namespace stan {
namespace math {

// Frechet(y|alpha, sigma)     [y > 0;  alpha > 0;  sigma > 0]
// FIXME: document
template <bool propto, typename T_y, typename T_shape, typename T_scale,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_shape, T_scale>* = nullptr>
return_type_t<T_y, T_shape, T_scale> frechet_lpdf(const T_y& y,
                                                  const T_shape& alpha,
                                                  const T_scale& sigma) {
  using T_partials_return = partials_return_t<T_y, T_shape, T_scale>;
  using T_y_ref = ref_type_t<T_y>;
  using T_alpha_ref = ref_type_t<T_shape>;
  using T_sigma_ref = ref_type_t<T_scale>;
  static const char* function = "frechet_lpdf";
  check_consistent_sizes(function, "Random variable", y, "Shape parameter",
                         alpha, "Scale parameter", sigma);
  T_y_ref y_ref = y;
  T_alpha_ref alpha_ref = alpha;
  T_sigma_ref sigma_ref = sigma;
  using std::pow;

  decltype(auto) y_val = to_ref(as_value_column_array_or_scalar(y_ref));
  decltype(auto) alpha_val = to_ref(as_value_column_array_or_scalar(alpha_ref));
  decltype(auto) sigma_val = to_ref(as_value_column_array_or_scalar(sigma_ref));

  check_positive(function, "Random variable", y_val);
  check_positive_finite(function, "Shape parameter", alpha_val);
  check_positive_finite(function, "Scale parameter", sigma_val);

  if (size_zero(y, alpha, sigma)) {
    return 0;
  }
  if (!include_summand<propto, T_y, T_shape, T_scale>::value) {
    return 0;
  }

  operands_and_partials<T_y_ref, T_alpha_ref, T_sigma_ref> ops_partials(
      y_ref, alpha_ref, sigma_ref);

  const auto& log_y
      = to_ref_if<(include_summand<propto, T_y, T_shape>::value
                   && include_summand<propto, T_shape, T_scale>::value)>(
          log(y_val));
  const auto& sigma_div_y_pow_alpha
      = to_ref_if<!is_constant_all<T_y, T_shape, T_scale>::value>(
          pow(sigma_val / y_val, alpha_val));

  size_t N = max_size(y, alpha, sigma);
  T_partials_return logp = -sum(sigma_div_y_pow_alpha);
  if (include_summand<propto, T_shape>::value) {
    logp += sum(log(alpha_val)) * N / math::size(alpha);
  }
  if (include_summand<propto, T_y, T_shape>::value) {
    logp -= sum((alpha_val + 1.0) * log_y) * N / max_size(y, alpha);
  }
  if (include_summand<propto, T_shape, T_scale>::value) {
    const auto& log_sigma
        = to_ref_if<!is_constant_all<T_shape>::value>(log(sigma_val));
    logp += sum(alpha_val * log_sigma) * N / max_size(alpha, sigma);
    if (!is_constant_all<T_shape>::value) {
      ops_partials.edge2_.partials_
          = inv(alpha_val) + (1 - sigma_div_y_pow_alpha) * (log_sigma - log_y);
    }
  }
  if (!is_constant_all<T_y>::value) {
    ops_partials.edge1_.partials_
        = (alpha_val * sigma_div_y_pow_alpha - (alpha_val + 1)) / y_val;
  }
  if (!is_constant_all<T_scale>::value) {
    ops_partials.edge3_.partials_
        = alpha_val / sigma_val * (1 - sigma_div_y_pow_alpha);
  }
  return ops_partials.build(logp);
}

template <typename T_y, typename T_shape, typename T_scale>
inline return_type_t<T_y, T_shape, T_scale> frechet_lpdf(const T_y& y,
                                                         const T_shape& alpha,
                                                         const T_scale& sigma) {
  return frechet_lpdf<false>(y, alpha, sigma);
}

}  // namespace math
}  // namespace stan
#endif
