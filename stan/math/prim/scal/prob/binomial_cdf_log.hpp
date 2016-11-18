#ifndef STAN_MATH_PRIM_SCAL_PROB_BINOMIAL_CDF_LOG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_BINOMIAL_CDF_LOG_HPP

#include <stan/math/prim/scal/prob/binomial_lcdf.hpp>

namespace stan {
  namespace math {

    /**
     * @deprecated use <code>binomial_lcdf</code>
     */
    template <typename T_n, typename T_N, typename T_prob>
    typename return_type<T_prob>::type
    binomial_cdf_log(const T_n& n, const T_N& N, const T_prob& theta) {
      return binomial_lcdf<T_n, T_N, T_prob>(n, N, theta);
    }

  }
}
#endif
