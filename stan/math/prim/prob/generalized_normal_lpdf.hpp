#ifndef STAN_MATH_PRIM_PROB_GENERALIZED_NORMAL_LPDF_HPP
#define STAN_MATH_PRIM_PROB_GENERALIZED_NORMAL_LPDF_HPP

#include <cmath>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <stan/math/prim/fun/abs.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/eval.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/inv_square.hpp>
#include <stan/math/prim/fun/lgamma.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/multiply_log.hpp>
#include <stan/math/prim/fun/pow.hpp>
#include <stan/math/prim/fun/promote_scalar.hpp>
#include <stan/math/prim/fun/select.hpp>
#include <stan/math/prim/fun/sign.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>
#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * The log of the generalized normal density for the specified scalar(s) given
 * the specified loc, scale and shape parameters. y, loc,
 * scale, or shape can each be either a scalar or a vector. Any vector inputs
 * must be the same length.
 *
 * <p>The result log probability is defined to be the sum of the
 * log probabilities for each y/mean/scale/shape tuple.
 *
 * @tparam T_y type of scalar
 * @tparam T_loc type of loc parameter
 * @tparam T_scale type of scale parameter
 * @tparam T_shape type of shape parameter
 * @param y (Sequence of) scalar(s)
 * @param loc (Sequence of) loc parameter(s)
 * @param scale (Sequence of) scale parameter(s)
 * @param shape (Sequence of) shape parameter(s)
 * @return The log of the product of the densities
 * @throw std::domain_error if scale or shape is not positive or not finite
 */
template <bool propto, typename T_y, typename T_loc, typename T_scale,
          typename T_shape,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_loc, T_scale, T_shape>* = nullptr>
inline return_type_t<T_y, T_loc, T_scale, T_shape> generalized_normal_lpdf(
    T_y&& y, T_loc&& loc, T_scale&& scale, T_shape&& shape) {
  using T_partials_return = partials_return_t<T_y, T_loc, T_scale, T_shape>;
  constexpr bool y_is_ad = is_autodiff_v<T_y>;
  constexpr bool loc_is_ad = is_autodiff_v<T_loc>;
  constexpr bool scale_is_ad = is_autodiff_v<T_scale>;
  constexpr bool shape_is_ad = is_autodiff_v<T_shape>;
  static constexpr const char* function = "generalized_normal_lpdf";
  check_consistent_sizes(function, "Random variable", y, "loc parameter", loc,
                         "Scale parameter", scale, "Shape parameter", shape);
  decltype(auto) y_ref = to_ref(std::forward<T_y>(y));
  decltype(auto) loc_ref = to_ref(std::forward<T_loc>(loc));
  decltype(auto) scale_ref = to_ref(std::forward<T_scale>(scale));
  decltype(auto) shape_ref = to_ref(std::forward<T_shape>(shape));
  auto y_val = as_value_column_array_or_scalar(y_ref);
  auto loc_val = as_value_column_array_or_scalar(loc_ref);
  auto scale_val = as_value_column_array_or_scalar(scale_ref);
  auto shape_val = as_value_column_array_or_scalar(shape_ref);

  check_not_nan(function, "Random variable", y_val);
  check_finite(function, "Location parameter", loc_val);
  check_positive_finite(function, "Scale parameter", scale_val);

  // With beta = +inf this could be defined to be uniform, but we don't support
  // that.
  check_positive_finite(function, "Shape parameter", shape_val);

  if (size_zero(y_ref, loc_ref, scale_ref, shape_ref)) {
    return 0;
  }
  if constexpr (propto && !is_any_autodiff_v<T_y, T_loc, T_scale, T_shape>) {
    return 0;
  } else {
    auto ops_partials
        = make_partials_propagator(y_ref, loc_ref, scale_ref, shape_ref);
    const size_t num_terms = max_size(y_ref, loc_ref, scale_ref, shape_ref);

    auto resid = eval((y_val - loc_val) / scale_val);
    auto resid_is_zero = as_array_or_scalar(value_of_rec(resid)) == 0;
    auto shape_eq_two = as_array_or_scalar(value_of_rec(shape_val)) == 2.0;
    // Replace zero bases with one so pow and multiply_log have finite tangents.
    auto residual_base = select(resid_is_zero, 1.0, abs(resid));
    auto generic_kernel = as_array_or_scalar(pow(residual_base, shape_val));
    auto normal_kernel
        = as_array_or_scalar(promote_scalar<T_partials_return>(square(resid)));

    // At shape == 2 the square expression preserves the normal Hessian at
    // the tie, whereas the generic expression would inherit abs'(0) == 0.
    auto kernel = eval(
        select(resid_is_zero,
               select(shape_eq_two, normal_kernel, T_partials_return(0.0)),
               generic_kernel));
    const T_partials_return kernel_log_density = -sum(kernel);

    if constexpr (y_is_ad || loc_is_ad || scale_is_ad) {
      auto inverse_scale = inv(scale_val);
      if constexpr (y_is_ad || loc_is_ad) {
        auto residual_sign = sign(resid);
        auto residual_power = pow(residual_base, shape_val - 1);
        auto generic_loc_derivative
            = residual_sign * shape_val * residual_power * inverse_scale;
        auto normal_loc_derivative = as_array_or_scalar(
            promote_scalar<T_partials_return>(2.0 * resid / scale_val));
        auto loc_derivative = select(
            resid_is_zero,
            select(shape_eq_two, normal_loc_derivative, T_partials_return(0.0)),
            generic_loc_derivative);
        if constexpr (y_is_ad) {
          partials<0>(ops_partials) = -loc_derivative;
        }
        if constexpr (loc_is_ad) {
          partials<1>(ops_partials) = loc_derivative;
        }
      }
      if constexpr (scale_is_ad) {
        partials<2>(ops_partials) = (shape_val * kernel - 1) * inverse_scale;
      }
    }

    T_partials_return log_density = kernel_log_density;
    if constexpr (!propto) {
      log_density -= LOG_TWO * num_terms;
    }
    if constexpr (!propto || scale_is_ad) {
      log_density -= sum(log(scale_val)) * (num_terms / math::size(scale));
    }
    if constexpr (!propto || shape_is_ad) {
      auto one_plus_inv_shape = inv(shape_val) + 1;
      if constexpr (shape_is_ad) {
        // The guarded base avoids NaN tangents from multiply_log at (0, 0).
        partials<3>(ops_partials)
            = digamma(one_plus_inv_shape) * inv_square(shape_val)
              - multiply_log(kernel, residual_base);
      }
      log_density
          -= sum(lgamma(one_plus_inv_shape)) * (num_terms / math::size(shape));
    }
    return ops_partials.build(log_density);
  }
}

template <typename T_y, typename T_loc, typename T_scale, typename T_shape>
inline return_type_t<T_y, T_loc, T_scale, T_shape> generalized_normal_lpdf(
    T_y&& y, T_loc&& loc, T_scale&& scale, T_shape&& shape) {
  return generalized_normal_lpdf<false>(
      std::forward<T_y>(y), std::forward<T_loc>(loc),
      std::forward<T_scale>(scale), std::forward<T_shape>(shape));
}

}  // namespace math
}  // namespace stan
#endif
