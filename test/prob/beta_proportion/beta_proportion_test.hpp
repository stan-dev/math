// Arguments: Doubles, Doubles, Doubles
#include <stan/math/prim/scal.hpp>

using stan::math::var;
using std::numeric_limits;
using std::vector;

class AgradDistributionsBetaProportion : public AgradDistributionTest {
 public:
  void valid_values(vector<vector<double> >& parameters,
                    vector<double>& log_prob) {
    vector<double> param(3);

    param[0] = 0.4;  // y
    param[1] = 0.5;  // mu (location)
    param[2] = 0.5;  // kappa (precision)
    parameters.push_back(param);
    log_prob.push_back(-0.9333428397413457);  // expected log_prob

    param[0] = 0.5;  // y
    param[1] = 0.2;  // mu (location)
    param[2] = 0.4;  // kappa (precision)
    parameters.push_back(param);
    log_prob.push_back(-1.6070080920051264);  // expected log_prob

    param[0] = 0.85;  // y
    param[1] = 0.15;  // mu (location)
    param[2] = 4.5;   // kappa (precision)
    parameters.push_back(param);
    log_prob.push_back(-4.7214376176246775);  // expected log_prob
  }

  void invalid_values(vector<size_t>& index, vector<double>& value) {
    // y
    index.push_back(0U);
    value.push_back(-1.0);

    index.push_back(0U);
    value.push_back(2.0);

    // mu
    index.push_back(1U);
    value.push_back(0.0);

    index.push_back(1U);
    value.push_back(-1.0);

    index.push_back(1U);
    value.push_back(numeric_limits<double>::infinity());

    index.push_back(1U);
    value.push_back(-numeric_limits<double>::infinity());

    // kappa
    index.push_back(2U);
    value.push_back(0.0);

    index.push_back(2U);
    value.push_back(-1.0);

    index.push_back(2U);
    value.push_back(numeric_limits<double>::infinity());

    index.push_back(2U);
    value.push_back(-numeric_limits<double>::infinity());
  }

  template <typename T_y, typename T_loc, typename T_prec, typename T3,
            typename T4, typename T5>
  typename stan::return_type<T_y, T_loc, T_prec>::type log_prob(
      const T_y& y, const T_loc& mu, const T_prec& kappa, const T3&, const T4&,
      const T5&) {
    return stan::math::beta_proportion_lpdf(y, mu, kappa);
  }

  template <bool propto, typename T_y, typename T_loc, typename T_prec,
            typename T3, typename T4, typename T5>
  typename stan::return_type<T_y, T_loc, T_prec>::type log_prob(
      const T_y& y, const T_loc& mu, const T_prec& kappa, const T3&, const T4&,
      const T5&) {
    return stan::math::beta_proportion_lpdf<propto>(y, mu, kappa);
  }

  template <typename T_y, typename T_loc, typename T_prec, typename T3,
            typename T4, typename T5>
  typename stan::return_type<T_y, T_loc, T_prec, T3, T4, T5>::type
  log_prob_function(const T_y& y, const T_loc& mu, const T_prec& kappa,
                    const T3&, const T4&, const T5&) {
    using stan::math::log1m;
    using std::log;
    return (mu * kappa - 1.0) * log(y) + ((1.0 - mu) * kappa - 1.0) * log1m(y)
           + lgamma(kappa) - lgamma(mu * kappa) - lgamma((1.0 - mu) * kappa);
  }
};
