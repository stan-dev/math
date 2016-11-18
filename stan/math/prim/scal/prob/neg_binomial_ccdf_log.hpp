#ifndef STAN_MATH_PRIM_SCAL_PROB_NEG_BINOMIAL_CCDF_LOG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_NEG_BINOMIAL_CCDF_LOG_HPP

#include <stan/math/prim/scal/prob/neg_binomial_lccdf.hpp>

namespace stan {
  namespace math {

    /**
     * @deprecated use <code>neg_binomial_lccdf</code>
     */
    template <typename T_n, typename T_shape,
              typename T_inv_scale>
    typename return_type<T_shape, T_inv_scale>::type
    neg_binomial_ccdf_log(const T_n& n, const T_shape& alpha,
                          const T_inv_scale& beta) {
      return neg_binomial_lccdf<T_n, T_shape, T_inv_scale>(n, alpha, beta);
    }

  }
}
#endif
