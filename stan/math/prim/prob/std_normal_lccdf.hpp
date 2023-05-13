#ifndef STAN_MATH_PRIM_PROB_STD_NORMAL_LCCDF_HPP
#define STAN_MATH_PRIM_PROB_STD_NORMAL_LCCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/erf.hpp>
#include <stan/math/prim/fun/erfc.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <cmath>

namespace stan {
namespace math {

template <
    typename T_y,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y>* = nullptr>
inline return_type_t<T_y> std_normal_lccdf(const T_y& y) {
  using T_partials_return = partials_return_t<T_y>;
  using std::exp;
  using std::log;
  using T_y_ref = ref_type_t<T_y>;
  static const char* function = "std_normal_lccdf";
  T_y_ref y_ref = y;
  check_not_nan(function, "Random variable", y_ref);

  if (size_zero(y)) {
    return 0;
  }

  T_partials_return lccdf(0.0);
  auto ops_partials = make_partials_propagator(y_ref);

  scalar_seq_view<T_y_ref> y_vec(y_ref);
  size_t N = stan::math::size(y);

  for (size_t n = 0; n < N; n++) {
    const T_partials_return y_dbl = y_vec.val(n);
    const T_partials_return scaled_y = y_dbl * INV_SQRT_TWO;

    T_partials_return one_m_erf;
    if (y_dbl < -37.5) {
      one_m_erf = 2.0;
    } else if (y_dbl < -5.0) {
      one_m_erf = 2.0 - erfc(-scaled_y);
    } else if (y_dbl > 8.25) {
      one_m_erf = 0.0;
    } else {
      one_m_erf = 1.0 - erf(scaled_y);
    }

    lccdf += LOG_HALF + log(one_m_erf);

    if (!is_constant_all<T_y>::value) {
      const T_partials_return rep_deriv
          = y_dbl > 8.25
                ? INFTY
                : SQRT_TWO_OVER_SQRT_PI * exp(-scaled_y * scaled_y) / one_m_erf;
      partials<0>(ops_partials)[n] -= rep_deriv;
    }
  }

  return ops_partials.build(lccdf);
}

}  // namespace math
}  // namespace stan
#endif
