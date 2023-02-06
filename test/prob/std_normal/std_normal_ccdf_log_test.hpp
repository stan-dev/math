// Arguments: Doubles
#include <stan/math/prim/prob/std_normal_lccdf.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/erf.hpp>
#include <stan/math/prim/fun/constants.hpp>

using std::vector;

class AgradCcdfLogStdNormal : public AgradCcdfLogTest {
 public:
  void valid_values(vector<vector<double> >& parameters,
                    vector<double>& ccdf_log) {
    vector<double> param(1);

    param[0] = 0;  // y
    parameters.push_back(param);
    ccdf_log.push_back(std::log(0.5));  // expected ccdf_log

    param[0] = 1;  // y
    parameters.push_back(param);
    ccdf_log.push_back(std::log(1.0 - 0.8413447460685429257765));

    param[0] = -2;  // y
    parameters.push_back(param);
    ccdf_log.push_back(std::log(1.0 - 0.02275013194817921205471));

    param[0] = -3.5;  // y
    parameters.push_back(param);
    ccdf_log.push_back(std::log(1.0 - 0.00023262907903552503635));
  }

  void invalid_values(vector<size_t>& index, vector<double>& value) {
    // y
  }

  bool has_lower_bound() { return false; }

  bool has_upper_bound() { return false; }

  template <typename T_y, typename T1, typename T2, typename T3, typename T4,
            typename T5>
  stan::return_type_t<T_y> ccdf_log(const T_y& y, const T1, const T2, const T3&,
                                    const T4&, const T5&) {
    return stan::math::std_normal_lccdf(y);
  }

  template <typename T_y, typename T1, typename T2, typename T3, typename T4,
            typename T5>
  stan::return_type_t<T_y> ccdf_log_function(const T_y& y, const T1, const T2,
                                             const T3&, const T4&, const T5&) {
    using std::log;
    return log(0.5 - 0.5 * erf(y * stan::math::INV_SQRT_TWO));
  }
};
