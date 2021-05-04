#ifndef STAN_MATH_PRIM_PROB_VON_MISES_LPDF_HPP
#define STAN_MATH_PRIM_PROB_VON_MISES_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/cos.hpp>
#include <stan/math/prim/fun/floor.hpp>
#include <stan/math/prim/fun/log_modified_bessel_first_kind.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/modified_bessel_first_kind.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/operands_and_partials.hpp>
#include <cmath>

namespace stan {
namespace math {

template <bool propto, typename T_y, typename T_loc, typename T_scale>
return_type_t<T_y, T_loc, T_scale> von_mises_lpdf(T_y const& y, T_loc const& mu,
                                                  T_scale const& kappa) {
  using T_partials_return = partials_return_t<T_y, T_loc, T_scale>;
  using T_y_ref = ref_type_if_t<!is_constant<T_y>::value, T_y>;
  using T_mu_ref = ref_type_if_t<!is_constant<T_loc>::value, T_loc>;
  using T_kappa_ref = ref_type_if_t<!is_constant<T_scale>::value, T_scale>;
  static char const* const function = "von_mises_lpdf";
  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Scale parameter", kappa);
  T_y_ref y_ref = y;
  T_mu_ref mu_ref = mu;
  T_kappa_ref kappa_ref = kappa;

  decltype(auto) y_val = to_ref(as_value_column_array_or_scalar(y_ref));
  decltype(auto) mu_val = to_ref(as_value_column_array_or_scalar(mu_ref));
  decltype(auto) kappa_val = to_ref(as_value_column_array_or_scalar(kappa_ref));
  check_finite(function, "Random variable", y_val);
  check_finite(function, "Location parameter", mu_val);
  check_nonnegative(function, "Scale parameter", kappa_val);
  check_finite(function, "Scale parameter", kappa_val);

  if (size_zero(y, mu, kappa)) {
    return 0;
  }
  if (!include_summand<propto, T_y, T_loc, T_scale>::value) {
    return 0;
  }

  operands_and_partials<T_y_ref, T_mu_ref, T_kappa_ref> ops_partials(
      y_ref, mu_ref, kappa_ref);

  const auto& cos_mu_minus_y
      = to_ref_if<!is_constant_all<T_scale>::value>(cos(mu_val - y_val));

  size_t N = max_size(y, mu, kappa);
  T_partials_return logp = sum(kappa_val * cos_mu_minus_y);
  if (include_summand<propto>::value) {
    logp -= LOG_TWO_PI * N;
  }
  if (include_summand<propto, T_scale>::value) {
    logp -= sum(log_modified_bessel_first_kind(0, kappa_val)) * N / size(kappa);
  }

  if (!is_constant_all<T_y, T_loc>::value) {
    const auto& sin_diff = sin(y_val - mu_val);
    auto kappa_sin
        = to_ref_if<(!is_constant_all<T_y>::value
                     && !is_constant_all<T_loc>::value)>(kappa_val * sin_diff);
    if (!is_constant_all<T_y>::value) {
      ops_partials.edge1_.partials_ = -kappa_sin;
    }
    if (!is_constant_all<T_loc>::value) {
      ops_partials.edge2_.partials_ = std::move(kappa_sin);
    }
  }
  if (!is_constant_all<T_scale>::value) {
    ops_partials.edge3_.partials_
        = cos_mu_minus_y
          - modified_bessel_first_kind(-1, kappa_val)
                / modified_bessel_first_kind(0, kappa_val);
  }

  return ops_partials.build(logp);
}

template <typename T_y, typename T_loc, typename T_scale>
inline return_type_t<T_y, T_loc, T_scale> von_mises_lpdf(T_y const& y,
                                                         T_loc const& mu,
                                                         T_scale const& kappa) {
  return von_mises_lpdf<false>(y, mu, kappa);
}

}  // namespace math
}  // namespace stan
#endif
