#ifndef STAN_MATH_PRIM_SCAL_PROB_SKEW_NORMAL_CCDF_LOG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_SKEW_NORMAL_CCDF_LOG_HPP

#include <stan/math/prim/scal/prob/skew_normal_lccdf.hpp>

namespace stan {
  namespace math {

    /**
     * @deprecated use <code>skew_normal_lccdf</code>
     */
    template <typename T_y, typename T_loc, typename T_scale, typename T_shape>
    typename return_type<T_y, T_loc, T_scale, T_shape>::type
    skew_normal_ccdf_log(const T_y& y, const T_loc& mu, const T_scale& sigma,
                         const T_shape& alpha) {
      return skew_normal_lccdf<T_y, T_loc,
                               T_scale, T_shape>(y, mu, sigma, alpha);
    }

  }
}
#endif
