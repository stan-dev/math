// Arguments: Ints, Doubles
#include <stan/math/prim/prob/yule_simon_lpmf.hpp>
#include <stan/math/prim/fun/lbeta.hpp>
#include <stan/math/prim/fun/lgamma.hpp>

using std::numeric_limits;
using std::vector;

class AgradDistributionsYuleSimon : public AgradDistributionTest {
 public:
  void valid_values(vector<vector<double> >& parameters,
                    vector<double>& log_prob) {
    vector<double> param(2);

    param[0] = 5;     // n
    param[1] = 20.0;  // alpha
    parameters.push_back(param);
    log_prob.push_back(-9.494202658325099);  // expected log_prob

    param[0] = 10;   // n
    param[1] = 5.5;  // alpha
    parameters.push_back(param);
    log_prob.push_back(-9.108616882863778);  // expected log_prob
  }

  void invalid_values(vector<size_t>& index, vector<double>& value) {
    // n
    index.push_back(0U);
    value.push_back(-1);

    index.push_back(0U);
    value.push_back(0);

    // alpha
    index.push_back(1U);
    value.push_back(0.0);

    index.push_back(1U);
    value.push_back(-1.0);

    index.push_back(1U);
    value.push_back(std::numeric_limits<double>::infinity());
  }

  template <class T_n, class T_alpha, typename T2, typename T3, typename T4,
            typename T5>
  stan::return_type_t<T_alpha> log_prob(const T_n& n, const T_alpha& alpha,
                                        const T2&, const T3&, const T4&,
                                        const T5&) {
    return stan::math::yule_simon_lpmf(n, alpha);
  }

  template <bool propto, class T_n, class T_alpha, typename T2, typename T3,
            typename T4, typename T5>
  stan::return_type_t<T_alpha> log_prob(const T_n& n, const T_alpha& alpha,
                                        const T2&, const T3&, const T4&,
                                        const T5&) {
    return stan::math::yule_simon_lpmf<propto>(n, alpha);
  }

  template <class T_n, class T_alpha, typename T2, typename T3, typename T4,
            typename T5>
  stan::return_type_t<T_alpha> log_prob_function(const T_n& n,
                                                 const T_alpha& alpha,
                                                 const T2&, const T3&,
                                                 const T4&, const T5&) {
    using stan::math::lbeta;
    using std::log;
    return log(alpha) + lbeta(n, alpha + 1.0);
  }
};
