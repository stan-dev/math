// Arguments: Ints, Doubles, Doubles
#include <stan/prob/distributions/univariate/discrete/zero_inflated.hpp>
#include <stan/prob/distributions/univariate/discrete/binomial.hpp>

#include <stan/math/functions/multiply_log.hpp>
#include <stan/math/functions/log1m.hpp>
#include <stan/math/functions/binomial_coefficient_log.hpp>

using std::vector;
using std::numeric_limits;
using stan::agrad::var;

class AgradDistributionsZIPoissonLog : public AgradDistributionTest {
public:
  void valid_values(vector<vector<double> >& parameters,
                    vector<double>& log_prob) {
    vector<double> param(3);

    //Density in R:
    //require(VGAM); eta <- 2.0; zeta <- -1.0; dzipois(0, exp(eta), 1/(1+exp(-zeta)), log=TRUE)
    param[0] = 0;           // n
    param[1] = 2.0;          // eta
    param[2] = -1.0;          // zeta
    parameters.push_back(param);
    log_prob.push_back(-1.31158325581604); // expected log_prob

    param[0] = 0;          // n
    param[1] = 3.1;          // eta
    param[2] = 4.1;          // zeta
    parameters.push_back(param);
    log_prob.push_back(-0.01643685); // expected log_prob
    
    
    param[0] = 0;          // n
    param[1] = 3.1;          // eta
    param[2] = -4.1;          // zeta
    parameters.push_back(param);
    log_prob.push_back(-4.116437); // expected log_prob

    param[0] = 1;          // n
    param[1] = 3.1;          // eta
    param[2] = -1.1;          // zeta
    parameters.push_back(param);
    log_prob.push_back(-19.38529); // expected log_prob
    
  }
 
  void invalid_values(vector<size_t>& index, 
                      vector<double>& value) {
    // n
    index.push_back(0U);
    value.push_back(-1);
    
    // eta
    //index.push_back(1U);
    //value.push_back(0);
    
    // zeta
    //index.push_back(2U);
    //value.push_back(0);
  }

  template <class T_n, class T_log_rate, class T_zi,
            typename T3, typename T4, typename T5, 
            typename T6, typename T7, typename T8, 
            typename T9>
  typename stan::return_type<T_log_rate, T_zi>::type 
  log_prob(const T_n& n, const T_log_rate& eta, const T_zi& zeta,
     const T3&, const T4&, const T5&, 
     const T6&, const T7&, const T8&, 
     const T9&) {
    return stan::prob::zi_poisson_log_log(n, eta, zeta);
  }

  template <bool propto, 
      class T_n, class T_log_rate, class T_zi,
            typename T3, typename T4, typename T5, 
            typename T6, typename T7, typename T8, 
            typename T9>
  typename stan::return_type<T_log_rate, T_zi>::type 
  log_prob(const T_n& n, const T_log_rate& eta, const T_zi& zeta,
     const T3&, const T4&, const T5&, 
     const T6&, const T7&, const T8&, 
     const T9&) {
    return stan::prob::zi_poisson_log_log<propto>(n, eta, zeta);
  }
  

  template <class T_n, class T_log_rate, class T_zi,
            typename T3, typename T4, typename T5, 
            typename T6, typename T7, typename T8, 
            typename T9>
  var log_prob_function(const T_n& n, const T_log_rate& eta, const T_zi& zeta,
      const T3&, const T4&, const T5&, 
      const T6&, const T7&, const T8&, 
      const T9&) {
    using stan::math::binomial_coefficient_log;
    using stan::math::log_sum_exp;
    using stan::math::log1p_exp;
    using boost::math::lgamma;
    using stan::prob::include_summand;
    using std::exp;
    using std::log;

    var logp(0);
    if (n==0) {
      if (include_summand<true, T_log_rate, T_zi>::value)
        logp += log_sum_exp(-exp(eta), zeta);
      if (include_summand<true, T_zi>::value)
        logp -= log1p_exp(zeta);
    } else {
      //logp += lgamma(n+1);
      if (include_summand<true, T_log_rate>::value) {
        logp -= exp(eta);
        logp += n*eta;      
      }
      if (include_summand<true, T_zi>::value)
        logp -= log1p_exp(zeta);
    }
    return logp;
  }
};

