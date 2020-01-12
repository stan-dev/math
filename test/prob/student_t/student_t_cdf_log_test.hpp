// Arguments: Doubles, Doubles, Doubles, Doubles
#include <stan/math/prim.hpp>

using stan::math::var;
using std::numeric_limits;
using std::vector;

class AgradCdfLogStudentT : public AgradCdfLogTest {
 public:
  void valid_values(vector<vector<double> >& parameters,
                    vector<double>& cdf_log) {
    vector<double> param(4);

    param[0] = 5.0;  // y
    param[1] = 1.5;  // nu (Degrees of Freedom)
    param[2] = 3.3;  // mu (Location)
    param[3] = 1.0;  // sigma (Scale)
    parameters.push_back(param);
    cdf_log.push_back(std::log(0.8646688779244795508561));  // expected CDF_log

    param[0] = 2.5;  // y
    param[1] = 3.5;  // nu (Degrees of Freedom)
    param[2] = 3.3;  // mu (Location)
    param[3] = 1.0;  // sigma (Scale)
    parameters.push_back(param);
    cdf_log.push_back(std::log(0.2372327883473262233327));  // expected CDF_log
  }

  void invalid_values(vector<size_t>& index, vector<double>& value) {
    // nu
    index.push_back(1U);
    value.push_back(-1.0);

    index.push_back(1U);
    value.push_back(0.0);

    index.push_back(1U);
    value.push_back(numeric_limits<double>::infinity());

    // mu
    index.push_back(2U);
    value.push_back(-numeric_limits<double>::infinity());

    index.push_back(2U);
    value.push_back(numeric_limits<double>::infinity());

    // sigma
    index.push_back(3U);
    value.push_back(-1.0);

    index.push_back(3U);
    value.push_back(0.0);

    index.push_back(3U);
    value.push_back(numeric_limits<double>::infinity());
  }

  bool has_lower_bound() { return false; }

  bool has_upper_bound() { return false; }

  template <typename T_y, typename T_dof, typename T_loc, typename T_scale,
            typename T4, typename T5>
  typename stan::return_type<T_y, T_dof, T_loc, T_scale>::type cdf_log(
      const T_y& y, const T_dof& nu, const T_loc& mu, const T_scale& sigma,
      const T4&, const T5&) {
    return stan::math::student_t_cdf_log(y, nu, mu, sigma);
  }

  template <typename T_y, typename T_dof, typename T_loc, typename T_scale,
            typename T4, typename T5>
  typename stan::return_type<T_y, T_dof, T_loc, T_scale>::type cdf_log_function(
      const T_y& y, const T_dof& nu, const T_loc& mu, const T_scale& sigma,
      const T4&, const T5&) {
    return stan::math::student_t_cdf_log(y, nu, mu, sigma);
  }
};
