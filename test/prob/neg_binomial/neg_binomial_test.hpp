// Arguments: Ints, Doubles, Doubles
#include <stan/math/prim.hpp>
#include <stdexcept>

using stan::math::var;
using std::numeric_limits;
using std::vector;

class AgradDistributionsNegBinomial : public AgradDistributionTest {
 public:
  void valid_values(vector<vector<double> >& parameters,
                    vector<double>& log_prob) {
    vector<double> param(3);

    param[0] = 10;   // n
    param[1] = 2.0;  // alpha
    param[2] = 1.5;  // beta
    parameters.push_back(param);
    log_prob.push_back(-7.786663293475162284096);  // expected log_prob

    param[0] = 100;  // n
    param[1] = 3.0;  // alpha
    param[2] = 3.5;  // beta
    parameters.push_back(param);
    log_prob.push_back(-142.6147368129045105434);  // expected log_prob

    param[0] = 10;
    param[1] = 1e6;
    param[2] = 1e5;
    parameters.push_back(param);
    log_prob.push_back(-2.0785666431081630812);  // expected log_prob
  }

  void invalid_values(vector<size_t>& index, vector<double>& value) {
    // n
    index.push_back(0U);
    value.push_back(-1);

    // alpha
    index.push_back(1U);
    value.push_back(0);

    // beta
    index.push_back(2U);
    value.push_back(0);
  }

  template <class T_n, class T_shape, class T_inv_scale, typename T3,
            typename T4, typename T5>
  stan::return_type_t<T_shape, T_inv_scale> log_prob(const T_n& n,
                                                     const T_shape& alpha,
                                                     const T_inv_scale& beta,
                                                     const T3&, const T4&,
                                                     const T5&) {
    return stan::math::neg_binomial_log(n, alpha, beta);
  }

  template <bool propto, class T_n, class T_shape, class T_inv_scale,
            typename T3, typename T4, typename T5>
  stan::return_type_t<T_shape, T_inv_scale> log_prob(const T_n& n,
                                                     const T_shape& alpha,
                                                     const T_inv_scale& beta,
                                                     const T3&, const T4&,
                                                     const T5) {
    return stan::math::neg_binomial_log<propto>(n, alpha, beta);
  }

  template <class T_n, class T_shape, class T_inv_scale, typename T3,
            typename T4, typename T5>
  stan::return_type_t<T_shape, T_inv_scale> log_prob_function(
      const T_n& n, const T_shape& alpha, const T_inv_scale& beta, const T3&,
      const T4&, const T5&) {
    using stan::math::binomial_coefficient_log;
    using stan::math::log1m;
    using stan::math::multiply_log;
    using std::log;

    if (n != 0)
      return binomial_coefficient_log<
                 typename stan::scalar_type<T_shape>::type>(n + alpha - 1.0, n)
             - n * log1p(beta) + alpha * log(beta / (1 + beta));
    return log(0.0);
  }
};
