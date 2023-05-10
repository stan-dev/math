// Arguments: Doubles, Doubles, Doubles, Doubles
#include <stan/math/prim/prob/skew_double_exponential_lccdf.hpp>
#include <stan/math/prim/fun/log1m_exp.hpp>
#include <stan/math/prim/fun/log1m.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/exp.hpp>

using stan::math::var;
using std::numeric_limits;
using std::vector;

class AgradCcdfLogSkewDoubleExponential : public AgradCcdfLogTest {
 public:
  void valid_values(vector<vector<double> >& parameters,
                    vector<double>& ccdf_log) {
    vector<double> param(4);

    param[0] = 0.0;  // y
    param[1] = 0.1;  // mu
    param[2] = 1.0;  // sigma
    param[3] = 0.5;  // skewness
    parameters.push_back(param);
    ccdf_log.push_back(-0.6022443516335636015668);  // expected log_ccdf

    param[0] = 1.0;  // y
    param[1] = 0.0;  // mu
    param[2] = 1.0;  // sigma
    param[3] = 0.5;  // skewness
    parameters.push_back(param);
    ccdf_log.push_back(-1.693147180559945175204);  // expected log_ccdf

    param[0] = -2.0;  // y
    param[1] = 0.0;   // mu
    param[2] = 1.0;   // sigma
    param[3] = 0.5;   // skewness
    parameters.push_back(param);
    ccdf_log.push_back(-0.07006592016028136138406);  // expected log_ccdf

    param[0] = -3.5;  // y
    param[1] = 1.9;   // mu
    param[2] = 7.2;   // sigma
    param[3] = 0.9;   // skewness
    parameters.push_back(param);
    ccdf_log.push_back(-1.490043636778108249175);  // expected log_ccdf

    param[0] = -3.5;  // y
    param[1] = 1.9;   // mu
    param[2] = 7.2;   // sigma
    param[3] = 0.1;   // skewness
    parameters.push_back(param);
    ccdf_log.push_back(-0.02626597639452391344528);  // expected log_ccdf

    param[0] = 3.5;  // y
    param[1] = 1.9;  // mu
    param[2] = 7.2;  // sigma
    param[3] = 0.9;  // skewness
    parameters.push_back(param);
    ccdf_log.push_back(-2.702585092994045368187);  // expected log_ccdf

    param[0] = 3.5;  // y
    param[1] = 1.9;  // mu
    param[2] = 7.2;  // sigma
    param[3] = 0.1;  // skewness
    parameters.push_back(param);
    ccdf_log.push_back(-0.1498049601022706511788);  // expected log_ccdf
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

  bool has_lower_bound() { return false; }

  bool has_upper_bound() { return false; }

  template <typename T_y, typename T_loc, typename T_scale, typename T_skewness,
            typename T4, typename T5>
  stan::return_type_t<T_y, T_loc, T_scale, T_skewness> ccdf_log(
      const T_y& y, const T_loc& mu, const T_scale& sigma,
      const T_skewness& tau, const T4&, const T5&) {
    return stan::math::skew_double_exponential_lccdf(y, mu, sigma, tau);
  }

  template <typename T_y, typename T_loc, typename T_scale, typename T_skewness,
            typename T4, typename T5>
  stan::return_type_t<T_y, T_loc, T_scale, T_skewness> ccdf_log_function(
      const T_y& y, const T_loc& mu, const T_scale& sigma,
      const T_skewness& tau, const T4&, const T5&) {
    using stan::math::log1m;
    using stan::math::log1m_exp;
    using std::exp;
    using std::log;

    if (y < mu) {
      return log1m_exp(log(tau) - 2 / sigma * (1 - tau) * (mu - y));
    } else {
      return log1m_exp(log1m((1 - tau) * exp(-2 / sigma * tau * (y - mu))));
    }
  }
};
