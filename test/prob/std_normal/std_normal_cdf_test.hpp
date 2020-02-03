// Arguments: Doubles
#include <stan/math/prim.hpp>

using std::vector;

class AgradCdfStdNormal : public AgradCdfTest {
 public:
  void valid_values(vector<vector<double>>& parameters, vector<double>& cdf) {
    vector<double> param(1);

    param[0] = 0.1;  // y
    parameters.push_back(param);
    cdf.push_back(0.5398278372770289814654);  // expected cdf

    param[0] = 1;  // y
    parameters.push_back(param);
    cdf.push_back(0.8413447460685429485852);  // expected cdf

    param[0] = -2;  // y
    parameters.push_back(param);
    cdf.push_back(0.0227501319481792072003);  // expected cdf

    param[0] = -3.5;  // y
    parameters.push_back(param);
    cdf.push_back(0.00023262907903552503635);  // expected cdf
  }

  void invalid_values(vector<size_t>& index, vector<double>& value) {
    // y
  }

  bool has_lower_bound() { return false; }

  bool has_upper_bound() { return false; }

  template <typename T_y, typename T1, typename T2, typename T3, typename T4,
            typename T5>
  stan::return_type_t<T_y> cdf(const T_y& y, const T1&, const T2&, const T3&,
                               const T4&, const T5&) {
    return stan::math::std_normal_cdf(y);
  }

  template <typename T_y, typename T1, typename T2, typename T3, typename T4,
            typename T5>
  stan::return_type_t<T_y> cdf_function(const T_y& y, const T1&, const T2&,
                                        const T3&, const T4&, const T5&) {
    using stan::math::SQRT_TWO;
    return (0.5 + 0.5 * erf(y / SQRT_TWO));
  }
};
