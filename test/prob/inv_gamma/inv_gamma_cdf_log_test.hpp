// Arguments: Doubles, Doubles, Doubles
#include <stan/math/prim.hpp>

using stan::math::var;
using std::numeric_limits;
using std::vector;

class AgradCdfLogInvGamma : public AgradCdfLogTest {
 public:
  void valid_values(vector<vector<double> >& parameters,
                    vector<double>& cdf_log) {
    vector<double> param(3);

    param[0] = 3.0;  // y
    param[1] = 0.5;  // alpha (Shape)
    param[2] = 3.3;  // beta (Scale)
    parameters.push_back(param);
    cdf_log.push_back(
        std::log(0.138010737568659559638900550956));  // expected cdf_log
  }

  void invalid_values(vector<size_t>& index, vector<double>& value) {
    // y
    index.push_back(0U);
    value.push_back(-1.0);

    // alpha
    index.push_back(1U);
    value.push_back(-1.0);

    index.push_back(1U);
    value.push_back(0.0);

    index.push_back(1U);
    value.push_back(numeric_limits<double>::infinity());

    // beta
    index.push_back(2U);
    value.push_back(-1.0);

    index.push_back(2U);
    value.push_back(0.0);

    index.push_back(2U);
    value.push_back(numeric_limits<double>::infinity());
  }

  bool has_lower_bound() { return true; }

  double lower_bound() { return 0.0; }

  bool has_upper_bound() { return false; }

  template <typename T_y, typename T_shape, typename T_scale, typename T3,
            typename T4, typename T5>
  typename stan::return_type<T_y, T_shape, T_scale>::type cdf_log(
      const T_y& y, const T_shape& alpha, const T_scale& beta, const T3&,
      const T4&, const T5&) {
    return stan::math::inv_gamma_cdf_log(y, alpha, beta);
  }

  template <typename T_y, typename T_shape, typename T_scale, typename T3,
            typename T4, typename T5>
  typename stan::return_type<T_y, T_shape, T_scale>::type cdf_log_function(
      const T_y& y, const T_shape& alpha, const T_scale& beta, const T3&,
      const T4&, const T5&) {
    using stan::math::gamma_q;

    return log(gamma_q(alpha, beta / y));
  }
};
