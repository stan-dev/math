// Arguments: Doubles, Doubles, Doubles
#include <stan/math/prim.hpp>

using stan::math::var;
using std::numeric_limits;
using std::vector;

class AgradCdfLogWeibull : public AgradCdfLogTest {
 public:
  void valid_values(vector<vector<double> >& parameters,
                    vector<double>& cdf_log) {
    vector<double> param(3);

    param[0] = 2.0;  // y
    param[1] = 1.0;  // alpha
    param[2] = 1.0;  // sigma
    parameters.push_back(param);
    cdf_log.push_back(std::log(0.8646647167633872976822));  // expected cdf_log

    param[0] = 0.25;  // y
    param[1] = 2.9;   // alpha
    param[2] = 1.8;   // sigma
    parameters.push_back(param);
    cdf_log.push_back(std::log(0.00325857114906742380625));  // expected cdf_log

    param[0] = 3.9;   // y
    param[1] = 1.7;   // alpha
    param[2] = 0.25;  // sigma
    parameters.push_back(param);
    cdf_log.push_back(std::log(1.0));  // expected cdf_log
  }

  void invalid_values(vector<size_t>& index, vector<double>& value) {
    // y
    index.push_back(0U);
    value.push_back(-1.0);

    index.push_back(0U);
    value.push_back(-numeric_limits<double>::infinity());

    // alpha
    index.push_back(1U);
    value.push_back(numeric_limits<double>::infinity());

    index.push_back(1U);
    value.push_back(-numeric_limits<double>::infinity());

    // sigma
    index.push_back(2U);
    value.push_back(0.0);

    index.push_back(2U);
    value.push_back(-1.0);
  }

  bool has_lower_bound() { return true; }

  double lower_bound() { return 0.0; }

  bool has_upper_bound() { return false; }

  template <typename T_y, typename T_shape, typename T_scale, typename T3,
            typename T4, typename T5>
  typename stan::return_type<T_y, T_shape, T_scale>::type cdf_log(
      const T_y& y, const T_shape& alpha, const T_scale& sigma, const T3&,
      const T4&, const T5&) {
    return stan::math::weibull_cdf_log(y, alpha, sigma);
  }

  template <typename T_y, typename T_shape, typename T_scale, typename T3,
            typename T4, typename T5>
  typename stan::return_type<T_y, T_shape, T_scale>::type cdf_log_function(
      const T_y& y, const T_shape& alpha, const T_scale& sigma, const T3&,
      const T4&, const T5&) {
    using std::log;
    using std::pow;
    return log(1.0 - exp(-pow(y / sigma, alpha)));
  }
};
