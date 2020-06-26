#ifndef STAN_MATH_PRIM_PROB_FRECHET_LCDF_HPP
#define STAN_MATH_PRIM_PROB_FRECHET_LCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log1m.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/multiply_log.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <boost/random/weibull_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <stan/math/prim/functor/operands_and_partials.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T_y, typename T_shape, typename T_scale>
return_type_t<T_y, T_shape, T_scale> frechet_lcdf(const T_y& y,
                                                  const T_shape& alpha,
                                                  const T_scale& sigma) {
  using T_partials_return = partials_return_t<T_y, T_shape, T_scale>;
  using T_y_ref = ref_type_t<T_y>;
  using T_alpha_ref = ref_type_t<T_shape>;
  using T_sigma_ref = ref_type_t<T_scale>;
  static const char* function = "frechet_lcdf";
  T_y_ref y_ref = y;
  T_alpha_ref alpha_ref = alpha;
  T_sigma_ref sigma_ref = sigma;

  const auto& y_col = as_column_vector_or_scalar(y_ref);
  const auto& alpha_col = as_column_vector_or_scalar(alpha_ref);
  const auto& sigma_col = as_column_vector_or_scalar(sigma_ref);

  const auto& y_arr = as_array_or_scalar(y_col);
  const auto& alpha_arr = as_array_or_scalar(alpha_col);
  const auto& sigma_arr = as_array_or_scalar(sigma_col);

  ref_type_t<decltype(value_of(y_arr))> y_val = value_of(y_arr);
  ref_type_t<decltype(value_of(alpha_arr))> alpha_val = value_of(alpha_arr);
  ref_type_t<decltype(value_of(sigma_arr))> sigma_val = value_of(sigma_arr);

  check_positive(function, "Random variable", y_val);
  check_positive_finite(function, "Shape parameter", alpha_val);
  check_positive_finite(function, "Scale parameter", sigma_val);

  if (size_zero(y, alpha, sigma)) {
    return 0.0;
  }

  operands_and_partials<T_y_ref, T_alpha_ref, T_sigma_ref> ops_partials(
      y_ref, alpha_ref, sigma_ref);

  const auto& pow_n = to_ref_if<!is_constant_all<T_y, T_shape, T_scale>::value>(
      pow(sigma_val / y_val, alpha_val));
  T_partials_return cdf_log = -sum(pow_n);

  if (!is_constant_all<T_y>::value) {
    ops_partials.edge1_.partials_ = pow_n * alpha_val / y_val;
  }
  if (!is_constant_all<T_shape>::value) {
    ops_partials.edge2_.partials_ = pow_n * log(y_val / sigma_val);
  }
  if (!is_constant_all<T_scale>::value) {
    ops_partials.edge3_.partials_ = -pow_n * alpha_val / sigma_val;
  }
  return ops_partials.build(cdf_log);
}

}  // namespace math
}  // namespace stan
#endif
