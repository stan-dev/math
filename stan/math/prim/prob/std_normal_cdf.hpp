#ifndef STAN_MATH_PRIM_PROB_STD_NORMAL_CDF_HPP
#define STAN_MATH_PRIM_PROB_STD_NORMAL_CDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/erf.hpp>
#include <stan/math/prim/fun/erfc.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/operands_and_partials.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Calculates the standard normal cumulative distribution function
 * for the given variate.
 *
 * \f$\Phi(x) = \frac{1}{\sqrt{2 \pi}} \int_{-\inf}^x e^{-t^2/2} dt\f$.
 *
 * @tparam T_y type of y
 * @param y scalar variate
 * @return The standard normal cdf evaluated at the specified argument.
 */
template <
    typename T_y,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T_y>* = nullptr>
inline return_type_t<T_y> std_normal_cdf(const T_y& y) {
  using T_partials_return = partials_return_t<T_y>;
  using std::exp;
  using T_y_ref = ref_type_t<T_y>;
  static const char* function = "std_normal_cdf";
  T_y_ref y_ref = y;
  check_not_nan(function, "Random variable", y_ref);

  if (size_zero(y)) {
    return 1.0;
  }

  T_partials_return cdf(1.0);
  operands_and_partials<T_y_ref> ops_partials(y_ref);

  scalar_seq_view<T_y_ref> y_vec(y_ref);
  size_t N = stan::math::size(y);

  for (size_t n = 0; n < N; n++) {
    const T_partials_return y_dbl = y_vec.val(n);
    const T_partials_return scaled_y = y_dbl * INV_SQRT_TWO;
    T_partials_return cdf_n;
    if (y_dbl < -37.5) {
      cdf_n = 0.0;
    } else if (y_dbl < -5.0) {
      cdf_n = 0.5 * erfc(-scaled_y);
    } else if (y_dbl > 8.25) {
      cdf_n = 1;
    } else {
      cdf_n = 0.5 * (1.0 + erf(scaled_y));
    }

    cdf *= cdf_n;

    if (!is_constant_all<T_y>::value) {
      const T_partials_return rep_deriv
          = (y_dbl < -37.5)
                ? 0.0
                : INV_SQRT_TWO_PI * exp(-scaled_y * scaled_y) / cdf_n;
      if (!is_constant_all<T_y>::value) {
        ops_partials.edge1_.partials_[n] += rep_deriv;
      }
    }
  }

  if (!is_constant_all<T_y>::value) {
    for (size_t n = 0; n < N; ++n) {
      ops_partials.edge1_.partials_[n] *= cdf;
    }
  }

  return ops_partials.build(cdf);
}

}  // namespace math
}  // namespace stan
#endif
