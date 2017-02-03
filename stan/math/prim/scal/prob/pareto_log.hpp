#ifndef STAN_MATH_PRIM_SCAL_PROB_PARETO_LOG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_PARETO_LOG_HPP

#include <stan/math/prim/scal/prob/pareto_lpdf.hpp>

namespace stan {
  namespace math {

    /**
     * @deprecated use <code>pareto_lpdf</code>
     */
    template <bool propto,
              typename T_y, typename T_scale, typename T_shape>
    typename return_type<T_y, T_scale, T_shape>::type
    pareto_log(const T_y& y, const T_scale& y_min, const T_shape& alpha) {
      return pareto_lpdf<propto, T_y, T_scale, T_shape>(y, y_min, alpha);
    }

    /**
     * @deprecated use <code>pareto_lpdf</code>
     */
    template <typename T_y, typename T_scale, typename T_shape>
    inline
    typename return_type<T_y, T_scale, T_shape>::type
    pareto_log(const T_y& y, const T_scale& y_min, const T_shape& alpha) {
      return pareto_lpdf<T_y, T_scale, T_shape>(y, y_min, alpha);
    }

  }
}
#endif
