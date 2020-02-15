#ifndef STAN_MATH_PRIM_PROB_FRECHET_LCCDF_HPP
#define STAN_MATH_PRIM_PROB_FRECHET_LCCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log1m.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/multiply_log.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <boost/random/weibull_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T_y, typename T_shape, typename T_scale>
return_type_t<T_y, T_shape, T_scale> frechet_lccdf(const T_y& y,
                                                   const T_shape& alpha,
                                                   const T_scale& sigma) {
  using T_partials_return = partials_return_t<T_y, T_shape, T_scale>;

  static const char* function = "frechet_lccdf";

  if (size_zero(y, alpha, sigma)) {
    return 0.0;
  }

  T_partials_return ccdf_log(0.0);
  check_positive(function, "Random variable", y);
  check_positive_finite(function, "Shape parameter", alpha);
  check_positive_finite(function, "Scale parameter", sigma);

  operands_and_partials<T_y, T_shape, T_scale> ops_partials(y, alpha, sigma);

  using std::exp;
  using std::log;
  using std::pow;

  scalar_seq_view<T_y> y_vec(y);
  scalar_seq_view<T_scale> sigma_vec(sigma);
  scalar_seq_view<T_shape> alpha_vec(alpha);
  size_t N = max_size(y, sigma, alpha);

  for (size_t n = 0; n < N; n++) {
    const T_partials_return y_dbl = value_of(y_vec[n]);
    const T_partials_return sigma_dbl = value_of(sigma_vec[n]);
    const T_partials_return alpha_dbl = value_of(alpha_vec[n]);
    const T_partials_return pow_n = pow(sigma_dbl / y_dbl, alpha_dbl);
    const T_partials_return exp_n = exp(-pow_n);

    ccdf_log += log1m(exp_n);

    const T_partials_return rep_deriv = pow_n / (1.0 / exp_n - 1);
    if (!is_constant_all<T_y>::value) {
      ops_partials.edge1_.partials_[n] -= alpha_dbl / y_dbl * rep_deriv;
    }
    if (!is_constant_all<T_shape>::value) {
      ops_partials.edge2_.partials_[n] -= log(y_dbl / sigma_dbl) * rep_deriv;
    }
    if (!is_constant_all<T_scale>::value) {
      ops_partials.edge3_.partials_[n] += alpha_dbl / sigma_dbl * rep_deriv;
    }
  }
  return ops_partials.build(ccdf_log);
}

}  // namespace math
}  // namespace stan
#endif
