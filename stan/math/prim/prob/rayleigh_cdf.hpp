#ifndef STAN_MATH_PRIM_PROB_RAYLEIGH_CDF_HPP
#define STAN_MATH_PRIM_PROB_RAYLEIGH_CDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T_y, typename T_scale,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_scale>* = nullptr>
return_type_t<T_y, T_scale> rayleigh_cdf(const T_y& y, const T_scale& sigma) {
  using T_partials_return = partials_return_t<T_y, T_scale>;
  using T_y_ref = ref_type_if_t<!is_constant<T_y>::value, T_y>;
  using T_sigma_ref = ref_type_if_t<!is_constant<T_scale>::value, T_scale>;
  static const char* function = "rayleigh_cdf";
  check_consistent_sizes(function, "Random variable", y, "Scale parameter",
                         sigma);

  T_y_ref y_ref = y;
  T_sigma_ref sigma_ref = sigma;

  decltype(auto) y_val = to_ref(as_value_column_array_or_scalar(y_ref));
  decltype(auto) sigma_val = to_ref(as_value_column_array_or_scalar(sigma_ref));

  check_nonnegative(function, "Random variable", y_val);
  check_positive(function, "Scale parameter", sigma_val);

  if (size_zero(y, sigma)) {
    return 1.0;
  }

  auto ops_partials = make_partials_propagator(y_ref, sigma_ref);

  const auto& inv_sigma
      = to_ref_if<!is_constant_all<T_scale>::value>(inv(sigma_val));
  const auto& inv_sigma_square
      = to_ref_if<!is_constant_all<T_y, T_scale>::value>(square(inv_sigma));
  const auto& exp_val = to_ref_if<!is_constant_all<T_y, T_scale>::value>(
      exp(-0.5 * square(y_val) * inv_sigma_square));

  T_partials_return cdf = prod(1 - exp_val);

  if (!is_constant_all<T_y, T_scale>::value) {
    const auto& common_deriv = to_ref_if<(!is_constant_all<T_y>::value
                                          && !is_constant_all<T_scale>::value)>(
        y_val * inv_sigma_square * exp_val / (1.0 - exp_val) * cdf);
    if (!is_constant_all<T_y>::value) {
      partials<0>(ops_partials) = common_deriv;
    }
    if (!is_constant_all<T_scale>::value) {
      partials<1>(ops_partials) = -y_val * inv_sigma * common_deriv;
    }
  }

  return ops_partials.build(cdf);
}

}  // namespace math
}  // namespace stan
#endif
