// Arguments: Ints, Ints, Ints
#include <stan/math/prim.hpp>

using stan::math::var;
using std::vector;

class AgradCcdfLogDiscreteRange : public AgradCcdfLogTest {
 public:
  void valid_values(vector<vector<double>>& parameters,
                    vector<double>& ccdf_log) {
    vector<double> param(3);

    param[0] = 3;  // y
    param[1] = 1;  // lower
    param[2] = 5;  // upper
    parameters.push_back(param);
    ccdf_log.push_back(log(2.0 / 5));  // expected ccdf_log

    param[0] = 9;   // y
    param[1] = 5;   // lower
    param[2] = 15;  // upper
    parameters.push_back(param);
    ccdf_log.push_back(log(6.0 / 11));  // expected ccdf_log

    param[0] = 0;   // y
    param[1] = -4;  // lower
    param[2] = 5;   // upper
    parameters.push_back(param);
    ccdf_log.push_back(log(5.0 / 10));  // expected ccdf_log
  }

  void invalid_values(vector<size_t>& /*index*/, vector<double>& /*value*/) {
    // y

    // lower

    // upper
  }

  bool has_lower_bound() { return false; }

  bool has_upper_bound() { return false; }

  template <class T_y, class T_lower, class T_upper, typename T3, typename T4,
            typename T5>
  stan::return_type_t<T_y, T_lower, T_upper> ccdf_log(const T_y& y,
                                                      const T_lower& lower,
                                                      const T_upper& upper,
                                                      const T3&, const T4&,
                                                      const T5&) {
    return stan::math::discrete_range_lccdf(y, lower, upper);
  }

  template <class T_y, class T_lower, class T_upper, typename T3, typename T4,
            typename T5>
  stan::return_type_t<T_y, T_lower, T_upper> ccdf_log_function(
      const T_y& y, const T_lower& lower, const T_upper& upper, const T3&,
      const T4&, const T5&) {
    if (y < lower) {
      return 0.0;
    }

    if (y >= upper) {
      return stan::math::LOG_ZERO;
    }

    return log((upper - y) / (upper - lower + 1));
  }
};
