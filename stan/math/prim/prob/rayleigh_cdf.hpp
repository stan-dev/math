#ifndef STAN_MATH_PRIM_PROB_RAYLEIGH_CDF_HPP
#define STAN_MATH_PRIM_PROB_RAYLEIGH_CDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/operands_and_partials.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T_y, typename T_scale>
return_type_t<T_y, T_scale> rayleigh_cdf(const T_y& y, const T_scale& sigma) {
  using T_partials_return = partials_return_t<T_y, T_scale>;
  using T_y_ref = ref_type_if_t<!is_constant<T_y>::value, T_y>;
  using T_sigma_ref = ref_type_if_t<!is_constant<T_scale>::value, T_scale>;
  static const char* function = "rayleigh_cdf";
  check_consistent_sizes(function, "Random variable", y, "Scale parameter",
                         sigma);

  T_y_ref y_ref = y;
  T_sigma_ref sigma_ref = sigma;

  const auto& y_col = as_column_vector_or_scalar(y_ref);
  const auto& sigma_col = as_column_vector_or_scalar(sigma_ref);

  const auto& y_arr = as_array_or_scalar(y_col);
  const auto& sigma_arr = as_array_or_scalar(sigma_col);

  ref_type_t<decltype(value_of(y_arr))> y_val = value_of(y_arr);
  ref_type_t<decltype(value_of(sigma_arr))> sigma_val = value_of(sigma_arr);

  check_nonnegative(function, "Random variable", y_val);
  check_positive(function, "Scale parameter", sigma_val);

  if (size_zero(y, sigma)) {
    return 1.0;
  }

  operands_and_partials<T_y_ref, T_sigma_ref> ops_partials(y_ref, sigma_ref);

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
      ops_partials.edge1_.partials_ = common_deriv;
    }
    if (!is_constant_all<T_scale>::value) {
      ops_partials.edge2_.partials_ = -y_val * inv_sigma * common_deriv;
    }
  }

  return ops_partials.build(cdf);
}

}  // namespace math
}  // namespace stan
#endif
