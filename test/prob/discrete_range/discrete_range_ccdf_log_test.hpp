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
    ccdf_log.push_back(-0.916290731874154995573);  // expected ccdf_log

    param[0] = 9;   // y
    param[1] = 5;   // lower
    param[2] = 15;  // upper
    parameters.push_back(param);
    ccdf_log.push_back(-0.606135803570315601085);  // expected ccdf_log

    param[0] = 0;   // y
    param[1] = -4;  // lower
    param[2] = 5;   // upper
    parameters.push_back(param);
    ccdf_log.push_back(-0.693147180559945286227);  // expected ccdf_log
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
    if (y < lower || y > upper) {
      return stan::math::LOG_ZERO;
    }

    return log((upper - y) / (upper - lower + 1));
  }
};
