// Arguments: Ints, Ints, Ints
#include <stan/math/prim/prob/discrete_range_lcdf.hpp>
#include <stan/math/prim/fun/log.hpp>

using stan::math::var;
using std::numeric_limits;
using std::vector;

class AgradCdfLogDiscreteRange : public AgradCdfLogTest {
 public:
  void valid_values(vector<vector<double>>& parameters,
                    vector<double>& cdf_log) {
    vector<double> param(3);

    param[0] = 3;  // y
    param[1] = 1;  // lower
    param[2] = 5;  // upper
    parameters.push_back(param);
    cdf_log.push_back(log(3.0 / 5));  // expected cdf_log

    param[0] = 9;   // y
    param[1] = 5;   // lower
    param[2] = 15;  // upper
    parameters.push_back(param);
    cdf_log.push_back(log(5.0 / 11));  // expected cdf_log

    param[0] = 0;   // y
    param[1] = -4;  // lower
    param[2] = 5;   // upper
    parameters.push_back(param);
    cdf_log.push_back(log(5.0 / 10));  // expected cdf_log
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
  double cdf_log(const T_y& y, const T_lower& lower, const T_upper& upper,
                 const T3&, const T4&, const T5&) {
    return stan::math::discrete_range_lcdf(y, lower, upper);
  }

  template <class T_y, class T_lower, class T_upper, typename T3, typename T4,
            typename T5>
  double cdf_log_function(const T_y& y, const T_lower& lower,
                          const T_upper& upper, const T3&, const T4&,
                          const T5&) {
    if (y < lower) {
      return stan::math::LOG_ZERO;
    }

    if (y > upper) {
      return 0.0;
    }

    return log((y - lower + 1.0) / (upper - lower + 1.0));
  }
};
