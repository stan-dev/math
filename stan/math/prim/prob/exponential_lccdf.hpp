#ifndef STAN_MATH_PRIM_PROB_EXPONENTIAL_LCCDF_HPP
#define STAN_MATH_PRIM_PROB_EXPONENTIAL_LCCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>

namespace stan {
namespace math {

template <typename T_y, typename T_inv_scale,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_inv_scale>* = nullptr>
return_type_t<T_y, T_inv_scale> exponential_lccdf(const T_y& y,
                                                  const T_inv_scale& beta) {
  using T_partials_return = partials_return_t<T_y, T_inv_scale>;
  using T_partials_array = Eigen::Array<T_partials_return, Eigen::Dynamic, 1>;
  using T_y_ref = ref_type_if_not_constant_t<T_y>;
  using T_beta_ref = ref_type_if_not_constant_t<T_inv_scale>;
  static const char* function = "exponential_lccdf";
  T_y_ref y_ref = y;
  T_beta_ref beta_ref = beta;

  auto y_val = to_ref(as_value_column_array_or_scalar(y_ref));
  auto beta_val = to_ref(as_value_column_array_or_scalar(beta_ref));

  check_nonnegative(function, "Random variable", y_val);
  check_positive_finite(function, "Inverse scale parameter", beta_val);

  if (size_zero(y, beta)) {
    return 0;
  }

  auto ops_partials = make_partials_propagator(y_ref, beta_ref);

  T_partials_return ccdf_log = -sum(beta_val * y_val);

  if (!is_constant_all<T_y>::value) {
    using beta_val_scalar = scalar_type_t<decltype(beta_val)>;
    using beta_val_array = Eigen::Array<beta_val_scalar, Eigen::Dynamic, 1>;
    if (is_vector<T_y>::value && !is_vector<T_inv_scale>::value) {
      partials<0>(ops_partials) = T_partials_array::Constant(
          math::size(y), -forward_as<beta_val_scalar>(beta_val));
    } else if (is_vector<T_inv_scale>::value) {
      partials<0>(ops_partials) = -forward_as<beta_val_array>(beta_val);
    } else {
      partials<0>(ops_partials)[0] = -sum(beta_val);
    }
  }
  if (!is_constant_all<T_inv_scale>::value) {
    using y_val_scalar = scalar_type_t<decltype(y_val)>;
    using y_val_array = Eigen::Array<y_val_scalar, Eigen::Dynamic, 1>;
    if (is_vector<T_inv_scale>::value && !is_vector<T_y>::value) {
      partials<1>(ops_partials) = T_partials_array::Constant(
          math::size(beta), -forward_as<y_val_scalar>(y_val));
    } else if (is_vector<T_y>::value) {
      partials<1>(ops_partials) = -forward_as<y_val_array>(y_val);
    } else {
      partials<1>(ops_partials)[0] = -sum(y_val);
    }
  }
  return ops_partials.build(ccdf_log);
}

}  // namespace math
}  // namespace stan
#endif
