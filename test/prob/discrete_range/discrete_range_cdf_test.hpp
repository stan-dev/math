// Arguments: Ints, Ints, Ints
#include <stan/math/prim.hpp>

using stan::math::var;
using std::numeric_limits;
using std::vector;

class AgradCdfDiscreteRange : public AgradCdfTest {
 public:
  void valid_values(vector<vector<double>>& parameters, vector<double>& cdf) {
    vector<double> param(3);

    param[0] = 3;  // y
    param[1] = 1;  // lower
    param[2] = 5;  // upper
    parameters.push_back(param);
    cdf.push_back(3.0 / 5);  // expected cdf

    param[0] = 9;   // y
    param[1] = 5;   // lower
    param[2] = 15;  // upper
    parameters.push_back(param);
    cdf.push_back(5.0 / 11);  // expected cdf

    param[0] = 0;   // y
    param[1] = -4;  // lower
    param[2] = 5;   // upper
    parameters.push_back(param);
    cdf.push_back(5.0 / 10);  // expected cdf
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
  double cdf(const T_y& y, const T_lower& lower, const T_upper& upper,
             const T3&, const T4&, const T5&) {
    return stan::math::discrete_range_cdf(y, lower, upper);
  }

  template <class T_y, class T_lower, class T_upper, typename T3, typename T4,
            typename T5>
  double cdf_function(const T_y& y, const T_lower& lower, const T_upper& upper,
                      const T3&, const T4&, const T5&) {
    if (y < lower) {
      return 0;
    }

    if (y >= upper) {
      return 1;
    }

    return (y - lower + 1.0) / (upper - lower + 1.0);
  }
};
