// Arguments: Doubles, Doubles, Doubles
#include <stan/math/prim.hpp>

using stan::math::var;
using std::numeric_limits;
using std::vector;

class AgradCcdfLogGamma : public AgradCcdfLogTest {
 public:
  void valid_values(vector<vector<double> >& parameters,
                    vector<double>& ccdf_log) {
    vector<double> param(3);

    param[0] = 1.0;  // y
    param[1] = 2.0;  // alpha
    param[2] = 2.0;  // beta
    parameters.push_back(param);
    ccdf_log.push_back(
        std::log(1.0 - 0.5939941502901618930466));  // expected ccdf_log

    param[0] = 2.0;   // y
    param[1] = 0.25;  // alpha
    param[2] = 0.75;  // beta
    parameters.push_back(param);
    ccdf_log.push_back(
        std::log(1.0 - 0.9665835558410209582547));  // expected ccdf_log

    param[0] = 1.0;  // y
    param[1] = 1.0;  // alpha
    param[2] = 1.0;  // beta
    parameters.push_back(param);
    ccdf_log.push_back(
        std::log(1.0 - 0.6321205588285576659757));  // expected ccdf_log

    param[0] = 16.0;  // y
    param[1] = 3.0;   // alpha
    param[2] = 3.0;   // beta
    parameters.push_back(param);
    ccdf_log.push_back(std::log(1.7116220633718619e-18));  // expected ccdf_log
  }

  void invalid_values(vector<size_t>& index, vector<double>& value) {
    // y
    index.push_back(0U);
    value.push_back(-1.0);

    // alpha
    index.push_back(1U);
    value.push_back(0.0);

    index.push_back(1U);
    value.push_back(-1.0);

    index.push_back(1U);
    value.push_back(numeric_limits<double>::infinity());

    index.push_back(1U);
    value.push_back(-numeric_limits<double>::infinity());

    // beta
    index.push_back(2U);
    value.push_back(0.0);

    index.push_back(2U);
    value.push_back(-1.0);

    index.push_back(2U);
    value.push_back(numeric_limits<double>::infinity());

    index.push_back(2U);
    value.push_back(-numeric_limits<double>::infinity());
  }

  bool has_lower_bound() { return true; }

  double lower_bound() { return 0.0; }

  bool has_upper_bound() { return false; }

  template <typename T_y, typename T_shape, typename T_inv_scale, typename T3,
            typename T4, typename T5>
  typename stan::return_type<T_y, T_shape, T_inv_scale>::type ccdf_log(
      const T_y& y, const T_shape& alpha, const T_inv_scale& beta, const T3&,
      const T4&, const T5&) {
    return stan::math::gamma_ccdf_log(y, alpha, beta);
  }

  template <typename T_y, typename T_shape, typename T_inv_scale, typename T3,
            typename T4, typename T5>
  typename stan::return_type<T_y, T_shape, T_inv_scale>::type ccdf_log_function(
      const T_y& y, const T_shape& alpha, const T_inv_scale& beta, const T3&,
      const T4&, const T5&) {
    using boost::math::gamma_q;
    using stan::math::gamma_q;
    using std::log;

    return log(gamma_q(alpha, beta * y));
  }
};
