#ifndef STAN_MATH_PRIM_PROB_STD_NORMAL_CDF_HPP
#define STAN_MATH_PRIM_PROB_STD_NORMAL_CDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/prod.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/std_normal_lcdf_impl.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Calculates the standard normal cumulative distribution function
 * for the given variate.
 *
 * \f$\Phi(x) = \frac{1}{\sqrt{2 \pi}} \int_{-\inf}^x e^{-t^2/2} dt\f$.
 *
 * Each factor is `0.5 * erfc(-x / sqrt(2))`, which has no cancellation for
 * either sign; the partials are the cdf times `phi(x) / Phi(x)`.
 *
 * @tparam T_y type of y
 * @param y scalar variate
 * @return The standard normal cdf evaluated at the specified argument.
 */
template <
    typename T_y,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y>* = nullptr>
inline return_type_t<T_y> std_normal_cdf(T_y&& y) {
  using T_y_ref = ref_type_if_not_constant_t<T_y>;
  static constexpr const char* function = "std_normal_cdf";
  T_y_ref y_ref = std::forward<T_y>(y);
  // to_ref evaluates an expression argument once for the check and kernel.
  decltype(auto) y_val = to_ref(as_value_column_array_or_scalar(y_ref));
  check_not_nan(function, "Random variable", y_val);

  if (size_zero(y_ref)) {
    return 1.0;
  }
  if constexpr (is_autodiff_v<T_y>) {
    const auto [values, ratios]
        = internal::std_normal_cdf_value_ratio<true>(y_val);
    const auto cdf = prod(values);
    auto ops_partials = make_partials_propagator(y_ref);
    partials<0>(ops_partials) = cdf * ratios;
    return ops_partials.build(cdf);
  } else {
    return prod(internal::std_normal_cdf_value_ratio<false>(y_val));
  }
}

}  // namespace math
}  // namespace stan
#endif
