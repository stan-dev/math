#ifndef STAN_MATH_PRIM_SCAL_PROB_EXP_MOD_NORMAL_LOG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_EXP_MOD_NORMAL_LOG_HPP

#include <stan/math/prim/scal/prob/exp_mod_normal_lpdf.hpp>

namespace stan {
  namespace math {

    /**
     * @deprecated use <code>exp_mod_normal_lpdf</code>
     */
    template <bool propto,
              typename T_y, typename T_loc, typename T_scale,
              typename T_inv_scale>
    typename return_type<T_y, T_loc, T_scale, T_inv_scale>::type
    exp_mod_normal_log(const T_y& y, const T_loc& mu, const T_scale& sigma,
                       const T_inv_scale& lambda) {
      return exp_mod_normal_lpdf<propto, T_y, T_loc,
                                 T_scale, T_inv_scale>(y, mu, sigma, lambda);
    }

    /**
     * @deprecated use <code>exp_mod_normal_lpdf</code>
     */
    template <typename T_y, typename T_loc, typename T_scale,
              typename T_inv_scale>
    inline
    typename return_type<T_y, T_loc, T_scale, T_inv_scale>::type
    exp_mod_normal_log(const T_y& y, const T_loc& mu, const T_scale& sigma,
                       const T_inv_scale& lambda) {
      return exp_mod_normal_lpdf<T_y, T_loc,
                                 T_scale, T_inv_scale>(y, mu, sigma, lambda);
    }

  }
}
#endif

