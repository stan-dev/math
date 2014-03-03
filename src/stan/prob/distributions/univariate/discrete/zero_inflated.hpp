#ifndef __STAN__PROB__DISTRIBUTIONS__UNIVARIATE__DISCRETE__POISSON_HPP__
#define __STAN__PROB__DISTRIBUTIONS__UNIVARIATE__DISCRETE__POISSON_HPP__

#include <boost/random/poisson_distribution.hpp>
#include <boost/random/variate_generator.hpp>

#include <limits>
#include <boost/math/special_functions/fpclassify.hpp>

#include <stan/agrad/partials_vari.hpp>
#include <stan/math/error_handling.hpp>
#include <stan/math/constants.hpp>
#include <stan/math/functions/multiply_log.hpp>
#include <stan/math/functions/value_of.hpp>
#include <stan/math/functions/log_sum_exp.hpp>
#include <stan/math/functions/log1p_exp.hpp>
#include <stan/meta/traits.hpp>
#include <stan/prob/traits.hpp>
#include <stan/prob/constants.hpp>

namespace stan {

  namespace prob {

    //ZIPoisson(n|exp(eta), zeta)
    template <bool propto,
              typename T_n, typename T_log_rate>
    typename return_type<T_log_rate>::type
    zi_poisson_log_log(const T_n& n, const T_log_rate& eta, const T_zi& zeta) {

      static const char* function = "stan::prob::zi_poisson_log_log(%1%)";
      
      using boost::math::lgamma;
      using stan::math::check_not_nan;
      using stan::math::check_nonnegative;
      using stan::math::value_of;
      using stan::math::check_consistent_sizes;
      using stan::prob::include_summand;
      using stan::math::multiply_log;
      using stan::math::log_sum_exp;
      using stan::math::log1p_exp;

      using std::exp;
      using std::log;
      
      // check if any vectors are zero length
      if (!(stan::length(n) && stan::length(eta)
            && stan::length(zeta)))
        return 0.0;

      // set up return value accumulator
      double logp(0.0);

      // validate args
      if (!check_nonnegative(function, n, "Random variable", &logp))
        return logp;
      if (!check_not_nan(function, eta,
                         "Log rate parameter", &logp))
        return logp;
      if (!(check_consistent_sizes(function,
                                   n, eta, "Random variable",
                                   "Log rate parameter", &logp)))
        return logp;
      if (!check_not_nan(function, zeta,
                         "Zero probability parameter", &logp))
        return logp;
      if (!check_bounded(function, zeta, 0.0, 1.0,
                         "Zero probability parameter", &logp))
        return logp;
      
      // check if no variables are involved and prop-to
      if (!include_summand<propto, T_log_rate, T_zi>::value)
        return 0.0;

      // set up expression templates wrapping scalars/vecs into vector views
      VectorView<const T_n> n_vec(n);
      VectorView<const T_log_rate> eta_vec(eta);
      VectorView<const T_zi> zeta_vec(zeta);
      size_t size = max_size(n, eta, zeta);

      // FIXME: first loop size of eta_vec, second loop if-ed for size==1
      for (size_t i = 0; i < size; i++)
        if (std::numeric_limits<double>::infinity() == eta_vec[i])
          return LOG_ZERO;
      for (size_t i = 0; i < size; i++)
        if (-std::numeric_limits<double>::infinity() == eta_vec[i] 
            && n_vec[i] != 0)
          return LOG_ZERO;
      
      // return accumulator with gradients
      agrad::OperandsAndPartials<T_log_rate, T_zi> operands_and_partials(eta, zeta);

      DoubleVectorView<true, is_vector<T_log_rate>::value>
        eta__(length(eta));
      for (size_t i = 0, size=length(eta); i < size; i++)
          eta__[i] = value_of(eta_vec[i]);

      DoubleVectorView<true, is_vector<T_zi>::value>
        zeta__(length(zeta));
      for (size_t i = 0, size=length(zeta); i < size; i++)
          zeta__[i] = value_of(zeta_vec[i]);

      for (size_t i = 0; i < size; i++) {
        
        if (n_vec[i] == 0) {
          
          if (include_summand<propto, T_log_rate, T_zi>::value)
            logp += log_sum_exp(-exp(eta__), zeta__);
          if (include_summand<propto, T_zi>::value)
            logp -=  log1p_exp(zeta__);

          // gradients
          if (!is_constant_struct<T_log_rate>::value)
            operands_and_partials.d_x1[i] -= 1
            / exp(log_sum_exp(-eta__, exp(eta__) + zeta__ - eta__));
          if (!is_constant_struct<T_zi>::value)
            operands_and_partials.d_x2[i] += 1/(1+exp(zeta__))
            - 1/(1+exp(exp(eta__)+zeta__));

        } else {
          
          if (include_summand<propto>::value)
            logp += lgamma(n_vec[i]+1);
          if (include_summand<propto, T_log_rate>::value) {
            logp -= exp(eta__);
            logp += n_vec[i]*eta__;
          }
          if (include_summand<propto, T_zi>::value)
            logp -=  log1p_exp(zeta__);

          // gradients
          if (!is_constant_struct<T_log_rate>::value)
            operands_and_partials.d_x1[i] += n_vec[i] - exp(eta__);
          if (!is_constant_struct<T_zi>::value)
            operands_and_partials.d_x2[i] -= 1/(1+exp(-zeta__));

        }
      }
      return operands_and_partials.to_var(logp);
    }
    
    template <typename T_n,
              typename T_log_rate>
    inline
    typename return_type<T_log_rate>::type
    poisson_log_log(const T_n& n, const T_log_rate& eta) {
      return poisson_log_log<false>(n, eta);
    }
  
  }
}
#endif
