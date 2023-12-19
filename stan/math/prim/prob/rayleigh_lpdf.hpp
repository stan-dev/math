#ifndef STAN_MATH_PRIM_PROB_RAYLEIGH_LPDF_HPP
#define STAN_MATH_PRIM_PROB_RAYLEIGH_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <cmath>

namespace stan {
namespace math {

template <bool propto, typename T_y, typename T_scale,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_scale>* = nullptr>
return_type_t<T_y, T_scale> rayleigh_lpdf(const T_y& y, const T_scale& sigma) {
  using T_partials_return = partials_return_t<T_y, T_scale>;
  using T_y_ref = ref_type_if_not_constant_t<T_y>;
  using T_sigma_ref = ref_type_if_not_constant_t<T_scale>;
  static const char* function = "rayleigh_lpdf";
  check_consistent_sizes(function, "Random variable", y, "Scale parameter",
                         sigma);

  T_y_ref y_ref = y;
  T_sigma_ref sigma_ref = sigma;

  decltype(auto) y_val = to_ref(as_value_column_array_or_scalar(y_ref));
  decltype(auto) sigma_val = to_ref(as_value_column_array_or_scalar(sigma_ref));

  check_positive(function, "Scale parameter", sigma_val);
  check_positive(function, "Random variable", y_val);

  if (size_zero(y, sigma)) {
    return 0.0;
  }
  if (!include_summand<propto, T_y, T_scale>::value) {
    return 0.0;
  }

  auto ops_partials = make_partials_propagator(y_ref, sigma_ref);

  const auto& inv_sigma
      = to_ref_if<!is_constant_all<T_y, T_scale>::value>(inv(sigma_val));
  const auto& y_over_sigma
      = to_ref_if<!is_constant_all<T_y, T_scale>::value>(y_val * inv_sigma);

  size_t N = max_size(y, sigma);
  T_partials_return logp = -0.5 * sum(square(y_over_sigma));
  if (include_summand<propto, T_scale>::value) {
    logp -= 2.0 * sum(log(sigma_val)) * N / math::size(sigma);
  }
  if (include_summand<propto, T_y>::value) {
    logp += sum(log(y_val)) * N / math::size(y);
  }

  if (!is_constant_all<T_y, T_scale>::value) {
    const auto& scaled_diff = to_ref_if<(!is_constant_all<T_y>::value
                                         && !is_constant_all<T_scale>::value)>(
        inv_sigma * y_over_sigma);
    if (!is_constant_all<T_y>::value) {
      partials<0>(ops_partials) = inv(y_val) - scaled_diff;
    }
    if (!is_constant_all<T_scale>::value) {
      edge<1>(ops_partials).partials_
          = y_over_sigma * scaled_diff - 2.0 * inv_sigma;
    }
  }

  return ops_partials.build(logp);
}

template <typename T_y, typename T_scale>
inline return_type_t<T_y, T_scale> rayleigh_lpdf(const T_y& y,
                                                 const T_scale& sigma) {
  return rayleigh_lpdf<false>(y, sigma);
}

}  // namespace math
}  // namespace stan
#endif
