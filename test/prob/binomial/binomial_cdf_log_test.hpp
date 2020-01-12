// Arguments: Ints, Ints, Doubles
#include <stan/math/prim.hpp>
#include <boost/math/special_functions/binomial.hpp>

using stan::math::var;
using std::numeric_limits;
using std::vector;

class AgradCdfLogBinomial : public AgradCdfLogTest {
 public:
  void valid_values(vector<vector<double> >& parameters,
                    vector<double>& cdf_log) {
    vector<double> param(3);

    param[0] = 17;   // Successes
    param[1] = 45;   // Trials
    param[2] = 0.5;  // Probability
    parameters.push_back(param);
    cdf_log.push_back(std::log(0.06757822542283530020679));  // expected cdf_log
  }

  void invalid_values(vector<size_t>& index, vector<double>& value) {
    // N (Trials)
    index.push_back(1U);
    value.push_back(-1);

    // p (Probability
    index.push_back(2U);
    value.push_back(-1e-4);

    index.push_back(2U);
    value.push_back(1 + 1e-4);
  }

  bool has_lower_bound() { return false; }

  bool has_upper_bound() { return false; }

  template <typename T_n, typename T_N, typename T_prob, typename T3,
            typename T4, typename T5>
  typename stan::return_type<T_prob>::type cdf_log(const T_n& n, const T_N& N,
                                                   const T_prob& theta,
                                                   const T3&, const T4&,
                                                   const T5&) {
    return stan::math::binomial_cdf_log(n, N, theta);
  }

  template <typename T_n, typename T_N, typename T_prob, typename T3,
            typename T4, typename T5>
  typename stan::return_type<T_prob>::type cdf_log_function(
      const T_n& n, const T_N& N, const T_prob& theta, const T3&, const T4&,
      const T5&) {
    using boost::math::binomial_coefficient;
    using std::exp;
    using std::log;

    typename stan::return_type<T_prob>::type cdf(0);

    for (int i = 0; i <= n; i++) {
      cdf += binomial_coefficient<double>(N, i)
             * exp(i * log(theta) + (N - i) * log(1 - theta));
    }

    return log(cdf);
  }
};
