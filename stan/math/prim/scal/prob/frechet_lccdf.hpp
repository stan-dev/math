#ifndef STAN_MATH_PRIM_SCAL_PROB_FRECHET_LCCDF_HPP
#define STAN_MATH_PRIM_SCAL_PROB_FRECHET_LCCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <boost/random/weibull_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/err/check_nonnegative.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/scal/err/check_positive.hpp>
#include <stan/math/prim/scal/err/check_positive_finite.hpp>
#include <stan/math/prim/scal/fun/size_zero.hpp>
#include <stan/math/prim/scal/fun/log1m.hpp>
#include <stan/math/prim/scal/fun/multiply_log.hpp>
#include <stan/math/prim/scal/fun/value_of.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T_y, typename T_shape, typename T_scale>
inline auto frechet_lccdf(const T_y& y, const T_shape& alpha,
                          const T_scale& sigma) {
  using T_partials = partials_return_t<T_y, T_shape, T_scale>;
  T_partials ccdf_log(0.0);


  using boost::math::tools::promote_args;
  using std::exp;
  using std::log;

  static const char* function = "frechet_lccdf";
  check_positive(function, "Random variable", y);
  check_positive_finite(function, "Shape parameter", alpha);
  check_positive_finite(function, "Scale parameter", sigma);

  const scalar_seq_view<T_y> y_vec(y);
  const scalar_seq_view<T_scale> sigma_vec(sigma);
  const scalar_seq_view<T_shape> alpha_vec(alpha);
  const size_t N = max_size(y, sigma, alpha);
  operands_and_partials<T_y, T_shape, T_scale> ops_partials(y, alpha, sigma);
  if (size_zero(y, alpha, sigma)) {
    return ops_partials.build(ccdf_log);
  }

  for (size_t n = 0; n < N; n++) {
    const T_partials y_dbl = value_of(y_vec[n]);
    const T_partials sigma_dbl = value_of(sigma_vec[n]);
    const T_partials alpha_dbl = value_of(alpha_vec[n]);
    const T_partials pow_ = pow(sigma_dbl / y_dbl, alpha_dbl);
    const T_partials exp_ = exp(-pow_);

    ccdf_log += log1m(exp_);

    const T_partials rep_deriv_ = pow_ / (1.0 / exp_ - 1);
    if (!is_constant_all<T_y>::value) {
      ops_partials.edge1_.partials_[n] -= alpha_dbl / y_dbl * rep_deriv_;
    }
    if (!is_constant_all<T_shape>::value) {
      ops_partials.edge2_.partials_[n] -= log(y_dbl / sigma_dbl) * rep_deriv_;
    }
    if (!is_constant_all<T_scale>::value) {
      ops_partials.edge3_.partials_[n] += alpha_dbl / sigma_dbl * rep_deriv_;
    }
  }
  return ops_partials.build(ccdf_log);
}

}  // namespace math
}  // namespace stan
#endif
