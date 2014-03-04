#ifndef __STAN__PROB__DISTRIBUTIONS__UNIVARIATE__DISCRETE__ZI_HPP__
#define __STAN__PROB__DISTRIBUTIONS__UNIVARIATE__DISCRETE__ZI_HPP__

#include <limits>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/math/special_functions/digamma.hpp>

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

    //ZIPoisson(n|exp(eta), inv_logit(zeta))
    template <bool propto,
              typename T_n, typename T_log_rate, typename T_zi>
    typename return_type<T_log_rate, T_zi>::type
    zi_poisson_log_log(const T_n& n, const T_log_rate& eta, const T_zi& zeta) {

      static const char* function = "stan::prob::zi_poisson_log_log(%1%)";
      
      using boost::math::lgamma;
      using stan::math::check_not_nan;
      using stan::math::check_nonnegative;
      using stan::math::value_of;
      using stan::math::check_consistent_sizes;
      using stan::math::check_finite;
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
      if (!check_finite(function, eta,
                         "Log rate parameter", &logp))
        return logp;
        
      if (!check_not_nan(function, zeta,
                         "Zero probability parameter", &logp))
        return logp;
      if (!check_finite(function, zeta,
                         "Zero probability parameter", &logp))
        return logp;
        
      if (!(check_consistent_sizes(function, n, eta, zeta,
                                   "Log rate parameter",
                                   "Random variable",
                                   "Zero probability parameter",
                                   &logp)))
        return logp;
      
      // check if no variables are involved and prop-to
      if (!include_summand<propto, T_log_rate, T_zi>::value)
        return 0.0;

      // set up expression templates wrapping scalars/vecs into vector views
      VectorView<const T_n> n_vec(n);
      VectorView<const T_log_rate> eta_vec(eta);
      VectorView<const T_zi> zeta_vec(zeta);
      size_t size = max_size(n, eta, zeta);
      
      // return accumulator with gradients
      agrad::OperandsAndPartials<T_log_rate, T_zi> operands_and_partials(eta, zeta);

      DoubleVectorView<true, is_vector<T_log_rate>::value>
        eta__(length(eta));
      for (size_t i = 0, size = length(eta); i < size; i++)
        eta__[i] = value_of(eta_vec[i]);

      DoubleVectorView<true, is_vector<T_zi>::value>
        zeta__(length(zeta));
      for (size_t i = 0, size = length(zeta); i < size; i++)
        zeta__[i] = value_of(zeta_vec[i]);

      for (size_t i = 0; i < size; i++) {
        
        if (n_vec[i] == 0) {
          
          if (include_summand<propto, T_log_rate, T_zi>::value)
            logp += log_sum_exp(-exp(eta__[i]), zeta__[i]);
          if (include_summand<propto, T_zi>::value)
            logp -=  log1p_exp(zeta__[i]);

          // gradients
          if (!is_constant_struct<T_log_rate>::value)
            operands_and_partials.d_x1[i] -= 1
            / exp(log_sum_exp(-eta__[i], exp(eta__[i]) + zeta__[i] - eta__[i]));
          if (!is_constant_struct<T_zi>::value)
            operands_and_partials.d_x2[i] += 1/(1+exp(zeta__[i]))
            - 1/(1+exp(exp(eta__[i])+zeta__[i]));

        } else {
          
          if (include_summand<propto>::value)
            logp += lgamma(n_vec[i]+1);
          if (include_summand<propto, T_log_rate>::value) {
            logp -= exp(eta__[i]);
            logp += n_vec[i]*eta__[i];
          }
          if (include_summand<propto, T_zi>::value)
            logp -=  log1p_exp(zeta__[i]);

          // gradients
          if (!is_constant_struct<T_log_rate>::value)
            operands_and_partials.d_x1[i] += n_vec[i] - exp(eta__[i]);
          if (!is_constant_struct<T_zi>::value)
            operands_and_partials.d_x2[i] -= 1/(1+exp(-zeta__[i]));

        }
      }
      return operands_and_partials.to_var(logp);
    }
    
    template <typename T_n,
              typename T_log_rate, typename T_zi>
    inline
    typename return_type<T_log_rate, T_zi>::type
    zi_poisson_log_log(const T_n& n, const T_log_rate& eta, const T_zi& zeta) {
      return zi_poisson_log_log<false>(n, eta, zeta);
    }













    //ZINegBin(n|exp(eta), phi, inv_logit(zeta))
    template <bool propto,
              typename T_n, typename T_log_location,
              typename T_inv_scale, typename T_zi>
    typename return_type<T_log_location, T_inv_scale, T_zi>::type
    zi_neg_binomial_log_log(const T_n& n, const T_log_location& eta,
    const T_inv_scale& phi, const T_zi& zeta) {

      static const char* function = "stan::prob::zi_neg_binomial_log_log(%1%)";
      
      using boost::math::lgamma;
      using stan::math::check_not_nan;
      using stan::math::check_nonnegative;
      using stan::math::value_of;
      using stan::math::check_consistent_sizes;
      using stan::math::check_positive;
      using stan::math::check_finite;
      using stan::prob::include_summand;
      using stan::math::multiply_log;
      using stan::math::log_sum_exp;
      using stan::math::log1p_exp;
      using boost::math::digamma;

      using std::exp;
      using std::log;
      
      // check if any vectors are zero length
      if (!(stan::length(n) && stan::length(eta)
            && stan::length(phi) && stan::length(zeta)))
        return 0.0;

      // set up return value accumulator
      double logp(0.0);

      // validate args
      if (!check_nonnegative(function, n, "Random variable", &logp))
        return logp;

      if (!check_not_nan(function, eta,
                         "Log location parameter", &logp))
        return logp;
      if (!check_finite(function, eta,
                         "Log location parameter", &logp))
        return logp;
        
      if (!check_not_nan(function, zeta,
                         "Zero probability parameter", &logp))
        return logp;
      if (!check_finite(function, zeta,
                         "Zero probability parameter", &logp))
        return logp;
        
      if (!check_not_nan(function, phi,
                         "Dispersion parameter", &logp))
        return logp;        
      if (!check_finite(function, phi,
                         "Dispersion parameter", &logp))
        return logp;
      if (!check_positive(function, phi,
                         "Dispersion parameter", &logp))
        return logp;
        
      if (!(check_consistent_sizes(function,
                                   n, eta, phi, zeta,
                                   "Random variable",
                                   "Log location parameter",
                                   "Dispersion parameter",
                                   "Zero probability parameter", &logp)))
        return logp;
      
      // check if no variables are involved and prop-to
      if (!include_summand<propto, T_log_location,
                           T_inv_scale, T_zi>::value)
        return 0.0;

      // set up expression templates wrapping scalars/vecs into vector views
      VectorView<const T_n> n_vec(n);
      VectorView<const T_log_location> eta_vec(eta);
      VectorView<const T_inv_scale> phi_vec(phi);
      VectorView<const T_zi> zeta_vec(zeta);
      size_t size = max_size(n, eta, phi, zeta);
      
      // return accumulator with gradients
      agrad::OperandsAndPartials<T_log_location, T_inv_scale, T_zi> 
                                  operands_and_partials(eta, phi, zeta);

      size_t len_ep = max_size(eta, phi);
      size_t len_np = max_size(n, phi);

      DoubleVectorView<true, is_vector<T_log_location>::value>
        eta__(length(eta));
      for (size_t i = 0, size = length(eta); i < size; i++)
        eta__[i] = value_of(eta_vec[i]);

      DoubleVectorView<true, is_vector<T_inv_scale>::value>
        phi__(length(phi));
      for (size_t i = 0, size = length(phi); i < size; i++)
        phi__[i] = value_of(phi_vec[i]);

      DoubleVectorView<true, is_vector<T_zi>::value>
        zeta__(length(zeta));
      for (size_t i = 0, size = length(zeta); i < size; i++)
        zeta__[i] = value_of(zeta_vec[i]);

      DoubleVectorView<true,is_vector<T_inv_scale>::value>
        log_phi(length(phi));
      for (size_t i = 0, size=length(phi); i < size; ++i)
        log_phi[i] = log(phi__[i]);

      DoubleVectorView<true,(is_vector<T_log_location>::value
                             || is_vector<T_inv_scale>::value)>
        logsumexp_eta_logphi(len_ep);
      for (size_t i = 0; i < len_ep; ++i)
        if (n_vec[i] != 0)
          logsumexp_eta_logphi[i] = log_sum_exp(eta__[i], log_phi[i]);

      DoubleVectorView<true,(is_vector<T_n>::value
                             || is_vector<T_inv_scale>::value)>
        n_plus_phi(len_np);
      for (size_t i = 0; i < len_np; ++i)
        if (n_vec[i] != 0)
          n_plus_phi[i] = n_vec[i] + phi__[i];

      for (size_t i = 0; i < size; i++) {
        
        if (n_vec[i] == 0) {
          double expression0 = log_phi[i] -
            log_sum_exp(eta__[i], log_phi[i]);
          double expression1 = (1 + phi__[i]) * expression0;
          double expression2 =
            log_sum_exp(zeta__[i], phi__[i] * expression0);
            
          if (include_summand<propto, T_log_location, T_inv_scale,
                              T_zi>::value)
            logp += expression2;
          if (include_summand<propto, T_zi>::value)
            logp -= log1p_exp(zeta__[i]);

          // gradients          
          if (!is_constant_struct<T_log_location>::value)
            operands_and_partials.d_x1[i] -= exp(eta__[i]
            + expression1 - expression2);
          if (!is_constant_struct<T_inv_scale>::value)
            operands_and_partials.d_x2[i] += exp(expression1
              - log_phi[i] - expression2)
              * (exp(eta__[i]) + (exp(eta__[i]) + phi__[i])*expression0);
          if (!is_constant_struct<T_zi>::value)
            operands_and_partials.d_x3[i] += exp(zeta__[i]
            - log1p_exp(zeta__[i]) - expression2)
            * (1 + exp(zeta__[i]) - exp(expression2));

        } else {
          
          if (include_summand<propto>::value)
            logp -= lgamma(n_vec[i] + 1.0);
          if (include_summand<propto, T_inv_scale>::value)
            logp += multiply_log(phi__[i], phi__[i]) - lgamma(phi__[i]);
          if (include_summand<propto, T_log_location,T_inv_scale>::value)
            logp -= (n_plus_phi[i])*logsumexp_eta_logphi[i];
          if (include_summand<propto, T_log_location>::value)
            logp += n_vec[i]*eta__[i];
          if (include_summand<propto, T_inv_scale>::value)
            logp += lgamma(n_plus_phi[i]);
  
          if (include_summand<propto, T_zi>::value)
            logp -=  log1p_exp(zeta__[i]);

          // gradients

          if (!is_constant_struct<T_log_location>::value)
            operands_and_partials.d_x1[i]
              += n_vec[i] - n_plus_phi[i]
              / (phi__[i]/exp(eta__[i]) + 1.0);
          if (!is_constant_struct<T_inv_scale>::value)
            operands_and_partials.d_x2[i]
              += 1.0 - n_plus_phi[i]/(exp(eta__[i]) + phi__[i])
              + log_phi[i] - logsumexp_eta_logphi[i]
              - digamma(phi__[i]) + digamma(n_plus_phi[i]);
          if (!is_constant_struct<T_zi>::value)
            operands_and_partials.d_x3[i] -= 1/(1+exp(-zeta__[i]));

        }
      }
      return operands_and_partials.to_var(logp);
    }
    
    template <typename T_n, typename T_log_location,
              typename T_inv_scale, typename T_zi>
    inline
    typename return_type<T_log_location, T_inv_scale, T_zi>::type
    zi_neg_binomial_log_log(const T_n& n, const T_log_location& eta,
    const T_inv_scale& phi, const T_zi& zeta) {
      return zi_neg_binomial_log_log<false>(n, eta, phi, zeta);
    }

  }
}
#endif
