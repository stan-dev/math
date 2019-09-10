#ifndef STAN_MATH_PRIM_SCAL_PROB_EXPONENTIAL_LCCDF_HPP
#define STAN_MATH_PRIM_SCAL_PROB_EXPONENTIAL_LCCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/err/check_nonnegative.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/scal/err/check_positive_finite.hpp>
#include <stan/math/prim/scal/fun/size_zero.hpp>
#include <stan/math/prim/scal/fun/value_of.hpp>

namespace stan {
namespace math {

template <typename T_y, typename T_inv_scale>
return_type_t<T_y, T_inv_scale> exponential_lccdf(const T_y& y,
                                                  const T_inv_scale& beta) {
  using T_partials_return = partials_return_t<T_y, T_inv_scale>;

  static const char* function = "exponential_lccdf";

  T_partials_return ccdf_log(0.0);
  if (size_zero(y, beta)) {
    return ccdf_log;
  }

  check_not_nan(function, "Random variable", y);
  check_nonnegative(function, "Random variable", y);
  check_positive_finite(function, "Inverse scale parameter", beta);

  operands_and_partials<T_y, T_inv_scale> ops_partials(y, beta);

  scalar_seq_view<T_y> y_vec(y);
  scalar_seq_view<T_inv_scale> beta_vec(beta);
  size_t N = max_size(y, beta);
  for (size_t n = 0; n < N; n++) {
    const T_partials_return beta_dbl = value_of(beta_vec[n]);
    const T_partials_return y_dbl = value_of(y_vec[n]);
    ccdf_log += -beta_dbl * y_dbl;

    if (!is_constant_all<T_y>::value) {
      ops_partials.edge1_.partials_[n] -= beta_dbl;
    }
    if (!is_constant_all<T_inv_scale>::value) {
      ops_partials.edge2_.partials_[n] -= y_dbl;
    }
  }
  return ops_partials.build(ccdf_log);
}

}  // namespace math
}  // namespace stan
#endif
