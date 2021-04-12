// Arguments: Doubles, Doubles, Doubles
#include <stan/math/prim.hpp>

using stan::math::var;
using std::numeric_limits;
using std::vector;

class AgradDistributionsLoglogistic : public AgradDistributionTest {
 public:
  void valid_values(vector<vector<double> >& parameters,
                    vector<double>& log_prob) {
    vector<double> param(3);

    param[0] = 1.0;  // y
    param[1] = 2.0;  // mu
    param[2] = 3.0;  // sigma
    parameters.push_back(param);
    log_prob.push_back(-1.21639532432449293253);  // expected log_prob

    param[0] = 1.4;  // y
    param[1] = 1.6;   // mu
    param[2] = 2.4;   // sigma
    parameters.push_back(param);
    log_prob.push_back(-0.87286484103664152556);  // expected log_prob
  }

  void invalid_values(vector<size_t>& index, vector<double>& value) {
    // y

    // mu
    index.push_back(1U);
    value.push_back(numeric_limits<double>::infinity());

    index.push_back(1U);
    value.push_back(-numeric_limits<double>::infinity());

    // sigma
    index.push_back(2U);
    value.push_back(0.0);

    index.push_back(2U);
    value.push_back(-1.0);

    index.push_back(2U);
    value.push_back(-numeric_limits<double>::infinity());

    index.push_back(2U);
    value.push_back(numeric_limits<double>::infinity());
  }

  template <typename T_y, typename T_loc, typename T_scale, typename T3,
            typename T4, typename T5>
  stan::return_type_t<T_y, T_loc, T_scale> log_prob(const T_y& y,
                                                    const T_loc& mu,
                                                    const T_scale& sigma,
                                                    const T3&, const T4&,
                                                    const T5&) {
    return stan::math::loglogistic_lpdf(y, mu, sigma);
  }

  template <bool propto, typename T_y, typename T_loc, typename T_scale,
            typename T3, typename T4, typename T5>
  stan::return_type_t<T_y, T_loc, T_scale> log_prob(const T_y& y,
                                                    const T_loc& mu,
                                                    const T_scale& sigma,
                                                    const T3&, const T4&,
                                                    const T5&) {
    return stan::math::loglogistic_lpdf<propto>(y, mu, sigma);
  }

  template <typename T_y, typename T_loc, typename T_scale, typename T3,
            typename T4, typename T5>
  stan::return_type_t<T_y, T_loc, T_scale> log_prob_function(
      const T_y& y, const T_loc& mu, const T_scale& sigma, const T3&, const T4&,
      const T5&) {
    using stan::math::log1p;
    return log(sigma) - log(mu) + (sigma - 1) * (log(y) - log(mu)) -
      2.0 * log1p(pow((y / mu), sigma));
    // T_partials_return logp = sum((log(sigma_val) - log(mu_val) + (sigma_val - 1) *
    //   (log(y_val) - log(mu_val)) -
    //   2 * log1p(pow((y_val / mu_val), sigma_val)));
    // return -(y - mu) / sigma - log(sigma) - 2.0 * log1p(exp(-(y - mu) / sigma));
  }
};
