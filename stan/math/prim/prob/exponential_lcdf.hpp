#ifndef STAN_MATH_PRIM_PROB_EXPONENTIAL_LCDF_HPP
#define STAN_MATH_PRIM_PROB_EXPONENTIAL_LCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/operands_and_partials.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T_y, typename T_inv_scale,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_inv_scale>* = nullptr>
return_type_t<T_y, T_inv_scale> exponential_lcdf(const T_y& y,
                                                 const T_inv_scale& beta) {
  using T_partials_return = partials_return_t<T_y, T_inv_scale>;
  using T_y_ref = ref_type_if_t<!is_constant<T_y>::value, T_y>;
  using T_beta_ref
      = ref_type_if_t<!is_constant<T_inv_scale>::value, T_inv_scale>;
  static const char* function = "exponential_lcdf";
  T_y_ref y_ref = y;
  T_beta_ref beta_ref = beta;

  decltype(auto) y_val = to_ref(as_value_column_array_or_scalar(y_ref));
  decltype(auto) beta_val = to_ref(as_value_column_array_or_scalar(beta_ref));

  check_nonnegative(function, "Random variable", y_val);
  check_positive_finite(function, "Inverse scale parameter", beta_val);

  if (size_zero(y, beta)) {
    return 0;
  }

  operands_and_partials<T_y_ref, T_beta_ref> ops_partials(y_ref, beta_ref);
  const auto& exp_val = to_ref_if<!is_constant_all<T_y, T_inv_scale>::value>(
      exp(-beta_val * y_val));

  T_partials_return cdf_log = sum(log1m(exp_val));

  if (!is_constant_all<T_y, T_inv_scale>::value) {
    const auto& rep_deriv = to_ref_if<(
        !is_constant_all<T_y>::value || !is_constant_all<T_inv_scale>::value)>(
        -exp_val / (1.0 - exp_val));
    if (!is_constant_all<T_y>::value) {
      ops_partials.edge1_.partials_ = -rep_deriv * beta_val;
    }
    if (!is_constant_all<T_inv_scale>::value) {
      ops_partials.edge2_.partials_ = -rep_deriv * y_val;
    }
  }
  return ops_partials.build(cdf_log);
}

}  // namespace math
}  // namespace stan
#endif
