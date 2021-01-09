// Arguments: Doubles, Doubles, Doubles, Doubles
#include <stan/math/prim.hpp>

using stan::math::var;
using std::numeric_limits;
using std::vector;

class AgradDistributionSkewDoubleExponential : public AgradDistributionTest {
 public:
  void valid_values(vector<vector<double> >& parameters,
                    vector<double>& log_prob) {
    vector<double> param(4);

    param[0] = 0.0;  // y
    param[1] = 0.1;  // mu
    param[2] = 1.0;  // sigma
    param[3] = 0.5;  // skewness
    parameters.push_back(param);
    log_prob.push_back(-0.7931471805599452640223);  // expected log_prob

    param[0] = 1.0;  // y
    param[1] = 0.0;  // mu
    param[2] = 1.0;  // sigma
    param[3] = 0.5;  // skewness
    parameters.push_back(param);
    log_prob.push_back(-1.693147180559945397249);  // expected log_prob

    param[0] = -2.0;  // y
    param[1] = 0.0;   // mu
    param[2] = 1.0;   // sigma
    param[3] = 0.5;   // skewness
    parameters.push_back(param);
    log_prob.push_back(-2.693147180559945397249);  // expected log_prob

    param[0] = -3.5;  // y
    param[1] = 1.9;   // mu
    param[2] = 7.2;   // sigma
    param[3] = 0.9;   // skewness
    parameters.push_back(param);
    log_prob.push_back(-3.838879454113936606774);  // expected log_prob

    param[0] = -3.5;  // y
    param[1] = 1.9;   // mu
    param[2] = 7.2;   // sigma
    param[3] = 0.1;   // skewness
    parameters.push_back(param);
    log_prob.push_back(-5.038879454113935452142);  // expected log_prob

    param[0] = 3.5;  // y
    param[1] = 1.9;  // mu
    param[2] = 7.2;  // sigma
    param[3] = 0.9;  // skewness
    parameters.push_back(param);
    log_prob.push_back(-4.088879454113937050863);  // expected log_prob

    param[0] = 3.5;  // y
    param[1] = 1.9;  // mu
    param[2] = 7.2;  // sigma
    param[3] = 0.1;  // skewness
    parameters.push_back(param);
    log_prob.push_back(-3.733323898558380093959);  // expected log_prob
  }

  void invalid_values(vector<size_t>& index, vector<double>& value) {
    // y

    // mu
    index.push_back(1U);
    value.push_back(-numeric_limits<double>::infinity());

    index.push_back(1U);
    value.push_back(numeric_limits<double>::infinity());

    // sigma
    index.push_back(2U);
    value.push_back(-1.0);

    index.push_back(2U);
    value.push_back(0.0);

    // skewness
    index.push_back(3U);
    value.push_back(-0.001);

    index.push_back(3U);
    value.push_back(1.001);
  }

  template <typename T_y, typename T_loc, typename T_scale, typename T_skewness,
            typename T4, typename T5>
  stan::return_type_t<T_y, T_loc, T_scale, T_skewness> log_prob(
      const T_y& y, const T_loc& mu, const T_scale& sigma,
      const T_skewness& tau, const T4&, const T5&) {
    return stan::math::skew_double_exponential_lpdf(y, mu, sigma, tau);
  }

  template <bool propto, typename T_y, typename T_loc, typename T_scale,
            typename T_skewness, typename T4, typename T5>
  stan::return_type_t<T_y, T_loc, T_scale, T_skewness> log_prob(
      const T_y& y, const T_loc& mu, const T_scale& sigma,
      const T_skewness& tau, const T4&, const T5&) {
    return stan::math::skew_double_exponential_lpdf<propto>(y, mu, sigma, tau);
  }

  template <typename T_y, typename T_loc, typename T_scale, typename T_skewness,
            typename T4, typename T5>
  stan::return_type_t<T_y, T_loc, T_scale, T_skewness> log_prob_function(
      const T_y& y, const T_loc& mu, const T_scale& sigma,
      const T_skewness& tau, const T4&, const T5&) {
    using stan::math::log1m;
    return log(2.0) + log(tau) + log1m(tau) - log(sigma)
           - 2.0 * ((y < mu) ? (1.0 - tau) * (mu - y) : tau * (y - mu)) / sigma;
  }
};
