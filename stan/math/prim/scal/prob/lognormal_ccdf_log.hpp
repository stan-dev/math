#ifndef STAN_MATH_PRIM_SCAL_PROB_LOGNORMAL_CCDF_LOG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_LOGNORMAL_CCDF_LOG_HPP

#include <stan/math/prim/scal/prob/lognormal_lccdf.hpp>

namespace stan {
  namespace math {

    /**
     * @deprecated use <code>lognormal_lccdf</code>
     */
    template <typename T_y, typename T_loc, typename T_scale>
    typename return_type<T_y, T_loc, T_scale>::type
    lognormal_ccdf_log(const T_y& y, const T_loc& mu, const T_scale& sigma) {
      return lognormal_lccdf<T_y, T_loc, T_scale>(y, mu, sigma);
    }

  }
}
#endif
