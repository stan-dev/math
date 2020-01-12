// Arguments: Doubles, Doubles
#include <stan/math/prim.hpp>

using stan::math::var;
using std::numeric_limits;
using std::vector;

class AgradCcdfLogInvChiSquare : public AgradCcdfLogTest {
 public:
  void valid_values(vector<vector<double> >& parameters,
                    vector<double>& ccdf_log) {
    vector<double> param(2);

    param[0] = 3.0;  // y
    param[1] = 0.5;  // Degrees of freedom
    parameters.push_back(param);
    ccdf_log.push_back(
        std::log(1.0 - 0.317528038297796704186230));  // expected ccdf_log
  }

  void invalid_values(vector<size_t>& index, vector<double>& value) {
    // y
    index.push_back(0U);
    value.push_back(-1.0);

    // nu
    index.push_back(1U);
    value.push_back(-1.0);

    index.push_back(1U);
    value.push_back(0.0);

    index.push_back(1U);
    value.push_back(numeric_limits<double>::infinity());
  }

  bool has_lower_bound() { return true; }

  double lower_bound() { return 0.0; }

  bool has_upper_bound() { return false; }

  template <typename T_y, typename T_dof, typename T2, typename T3, typename T4,
            typename T5>
  typename stan::return_type<T_y, T_dof>::type ccdf_log(const T_y& y,
                                                        const T_dof& nu,
                                                        const T2&, const T3&,
                                                        const T4&, const T5&) {
    return stan::math::inv_chi_square_ccdf_log(y, nu);
  }

  template <typename T_y, typename T_dof, typename T2, typename T3, typename T4,
            typename T5>
  typename stan::return_type<T_y, T_dof>::type ccdf_log_function(
      const T_y& y, const T_dof& nu, const T2&, const T3&, const T4&,
      const T5&) {
    using stan::math::gamma_p;

    return log(1.0 - gamma_p(0.5 * nu, 0.5 / y));
  }
};
