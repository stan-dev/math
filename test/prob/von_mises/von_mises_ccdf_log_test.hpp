// Arguments: Doubles, Doubles, Doubles
#include <stan/math/prim.hpp>

using stan::math::var;
using std::numeric_limits;
using std::vector;

class AgradCcdfLogVonMises : public AgradCcdfLogTest {
 public:
  void valid_values(vector<vector<double> >& parameters,
                    vector<double>& log_ccdf) {
    vector<double> param(3);

    param[0] = -1.5707463267948965;     // y
    param[1] = 6.283185307179586;       // mu
    param[2] = 1.0;                     // kappa
    parameters.push_back(param);
    log_ccdf.push_back(std::log(1.0 - 0.1097601896880533));  // expected Log_CCDF
  }

  void invalid_values(vector<size_t>& index, vector<double>& value) {
    // y
    index.push_back(0U);
    value.push_back(-stan::math::pi() - 1.0);

    index.push_back(0U);
    value.push_back(stan::math::pi() + 1.0);

    // mu
    index.push_back(1U);
    value.push_back(-numeric_limits<double>::infinity());

    index.push_back(1U);
    value.push_back(numeric_limits<double>::infinity());

    // kappa
    index.push_back(2U);
    value.push_back(-1.0);

    index.push_back(2U);
    value.push_back(-numeric_limits<double>::infinity());
  }

  bool has_lower_bound() { return true; }

  double lower_bound() { return -stan::math::pi(); }

  bool has_upper_bound() { return true; }

  double upper_bound() { return stan::math::pi(); }

  template <typename T_y, typename T_mu, typename T_scale,
            typename T3, typename T4, typename T5>
  stan::return_type_t<T_y, T_mu, T_scale> ccdf_log(
      const T_y& y, const T_mu& mu, const T_scale& scale,
      const T3&, const T4&, const T5&) {
    return stan::math::von_mises_ccdf_log(y, mu, scale);
  }

  template <typename T_y, typename T_mu, typename T_scale,
            typename T3, typename T4, typename T5>
  stan::return_type_t<T_y, T_mu, T_scale> ccdf_log_function(
      const T_y& y, const T_mu& mu, const T_scale& scale,
      const T3&, const T4&, const T5&) {
    return stan::math::von_mises_ccdf_log(y, mu, scale);
  }
};
