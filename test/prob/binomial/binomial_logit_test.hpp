// Arguments: Ints, Ints, Doubles
#include <stan/math/prim/prob/binomial_logit_lpmf.hpp>

using stan::math::var;
using std::numeric_limits;
using std::vector;

class AgradDistributionsBinomialLogit : public AgradDistributionTest {
 public:
  void valid_values(vector<vector<double> >& parameters,
                    vector<double>& log_prob) {
    vector<double> param(3);
    using stan::math::logit;

    param[0] = 10;          // n
    param[1] = 20;          // N
    param[2] = logit(0.4);  // alpha
    parameters.push_back(param);
    log_prob.push_back(-2.144372241799002765106);  // expected log_prob

    param[0] = 5;           // n
    param[1] = 15;          // N
    param[2] = logit(0.8);  // alpha
    parameters.push_back(param);
    log_prob.push_back(-9.202729812928724939525);  // expected log_prob
  }

  void invalid_values(vector<size_t>& index, vector<double>& value) {
    // n
    index.push_back(0U);
    value.push_back(-1);

    // N
    index.push_back(1U);
    value.push_back(-1);

    // alpha
    index.push_back(2U);
    value.push_back(std::numeric_limits<double>::infinity());
  }

  template <class T_n, class T_N, class T_prob, typename T3, typename T4,
            typename T5>
  stan::return_type_t<T_prob> log_prob(const T_n& n, const T_N& N,
                                       const T_prob& alpha, const T3&,
                                       const T4&, const T5&) {
    return stan::math::binomial_logit_lpmf(n, N, alpha);
  }

  template <bool propto, class T_n, class T_N, class T_prob, typename T3,
            typename T4, typename T5>
  stan::return_type_t<T_prob> log_prob(const T_n& n, const T_N& N,
                                       const T_prob& alpha, const T3&,
                                       const T4&, const T5&) {
    return stan::math::binomial_logit_lpmf<propto>(n, N, alpha);
  }

  template <class T_n, class T_N, class T_prob, typename T3, typename T4,
            typename T5>
  stan::return_type_t<T_prob> log_prob_function(const T_n& n, const T_N& N,
                                                const T_prob& alpha, const T3&,
                                                const T4&, const T5&) {
    using stan::math::binomial_coefficient_log;
    using stan::math::inv_logit;
    using stan::math::log1m;
    using stan::math::multiply_log;
    using std::log;

    return binomial_coefficient_log(N, n) + multiply_log(n, inv_logit(alpha))
           + (N - n) * log1m(inv_logit(alpha));
  }
};
