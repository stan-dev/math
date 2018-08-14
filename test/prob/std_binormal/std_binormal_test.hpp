// Arguments: Doubles, Doubles, Doubles
#include <stan/math/prim/scal.hpp>

using stan::math::var;
using std::numeric_limits;
using std::vector;

class AgradCdfLogBinormal : public AgradCdfLogTest {
 public:
  void valid_values(vector<vector<double> >& parameters, vector<double>& log_cdf) {
    vector<double> param(3);

    param[0] = 0.3;  // y_1
    param[1] = 1.4;  // y_2 
    param[2] = 0.7;  // rho 
    parameters.push_back(param);
    log_cdf.push_back(-0.4938715763653356);  // expected log_prob

    param[0] = 1.3;  // y_1
    param[1] = 0.4;  // y_2 
    param[2] = 0.3;  // rho 
    parameters.push_back(param);
    log_cdf.push_back(-0.4906122253047989);  // expected log_prob

    param[0] = 2.3;  // y_1
    param[1] = 3.4;  // y_2 
    param[2] = 0.3;  // rho 
    parameters.push_back(param);
    log_cdf.push_back(-0.01108711762415265);  // expected log_prob

    param[0] = 2.3;  // y_1
    param[1] = -3.4;  // y_2 
    param[2] = -0.3;  // rho 
    parameters.push_back(param);
    log_cdf.push_back(-8.105838688329269);  // expected log_prob
  }

  void invalid_values(vector<size_t>& index, vector<double>& value) {
    // y

    // rho 
    index.push_back(2U);
    value.push_back(numeric_limits<double>::infinity());

    index.push_back(2U);
    value.push_back(-numeric_limits<double>::infinity());

    index.push_back(2U);
    value.push_back(-1.3);

    index.push_back(2U);
    value.push_back(-1);

    index.push_back(2U);
    value.push_back(1);

    index.push_back(2U);
    value.push_back(1.3);
  }

  int random_variable_dimension() { return 2; }

  template <typename T_y_1, typename T_y_2, typename T_rho, typename T3,
            typename T4, typename T5>
  typename stan::return_type<T_y_1, T_y_2, T_rho>::type cdf_log(
      const T_y_1& y_1, const T_y_2& y_2, const T_rho& rho, const T3&, const T4&,
      const T5&) {
    return stan::math::std_binormal_lcdf(y_1, y_2, rho);
  }

  template <typename T_y_1, typename T_y_2, typename T_rho, typename T3,
            typename T4, typename T5>
  typename stan::return_type<T_y_1, T_y_2, T_rho, T3, T4, T5>::type
  cdf_log_function(const T_y_1& y_1, const T_y_2& y_2, const T_rho& rho,
                    const T3&, const T4&, const T5&) {
    using stan::math::binormal_integral_owens;
    return binormal_integral_owens(y_1, y_2, rho);
  }
};
