#ifndef STAN_MATH_PRIM_PROB_FRECHET_CDF_HPP
#define STAN_MATH_PRIM_PROB_FRECHET_CDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log1m.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/multiply_log.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <boost/random/weibull_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T_y, typename T_shape, typename T_scale,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_shape, T_scale>* = nullptr>
return_type_t<T_y, T_shape, T_scale> frechet_cdf(const T_y& y,
                                                 const T_shape& alpha,
                                                 const T_scale& sigma) {
  using T_partials_return = partials_return_t<T_y, T_shape, T_scale>;
  using T_y_ref = ref_type_t<T_y>;
  using T_alpha_ref = ref_type_t<T_shape>;
  using T_sigma_ref = ref_type_t<T_scale>;
  using std::pow;
  static constexpr const char* function = "frechet_cdf";
  T_y_ref y_ref = y;
  T_alpha_ref alpha_ref = alpha;
  T_sigma_ref sigma_ref = sigma;

  check_positive(function, "Random variable", y_ref);
  check_positive_finite(function, "Shape parameter", alpha_ref);
  check_positive_finite(function, "Scale parameter", sigma_ref);

  if (size_zero(y_ref, alpha_ref, sigma_ref)) {
    return 1.0;
  }

  T_partials_return cdf(1.0);
  auto ops_partials = make_partials_propagator(y_ref, alpha_ref, sigma_ref);

  scalar_seq_view<T_y> y_vec(y_ref);
  scalar_seq_view<T_scale> sigma_vec(sigma_ref);
  scalar_seq_view<T_shape> alpha_vec(alpha_ref);
  size_t N = max_size(y_ref, sigma_ref, alpha_ref);

  for (size_t n = 0; n < N; n++) {
    const T_partials_return y_dbl = y_vec.val(n);
    const T_partials_return sigma_dbl = sigma_vec.val(n);
    const T_partials_return alpha_dbl = alpha_vec.val(n);
    const T_partials_return pow_n = pow(sigma_dbl / y_dbl, alpha_dbl);
    const T_partials_return cdf_n = exp(-pow_n);

    cdf *= cdf_n;

    if (!is_constant_all<T_y>::value) {
      partials<0>(ops_partials)[n] += pow_n * alpha_dbl / y_dbl;
    }
    if (!is_constant_all<T_shape>::value) {
      partials<1>(ops_partials)[n] += pow_n * log(y_dbl / sigma_dbl);
    }
    if (!is_constant_all<T_scale>::value) {
      partials<2>(ops_partials)[n] -= pow_n * alpha_dbl / sigma_dbl;
    }
  }

  if (!is_constant_all<T_y>::value) {
    for (size_t n = 0; n < stan::math::size(y); ++n) {
      partials<0>(ops_partials)[n] *= cdf;
    }
  }
  if (!is_constant_all<T_shape>::value) {
    for (size_t n = 0; n < stan::math::size(alpha); ++n) {
      partials<1>(ops_partials)[n] *= cdf;
    }
  }
  if (!is_constant_all<T_scale>::value) {
    for (size_t n = 0; n < stan::math::size(sigma); ++n) {
      partials<2>(ops_partials)[n] *= cdf;
    }
  }
  return ops_partials.build(cdf);
}

}  // namespace math
}  // namespace stan
#endif
