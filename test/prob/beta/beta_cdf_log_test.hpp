// Arguments: Doubles, Doubles, Doubles
#include <stan/math/prim/prob/beta_cdf_log.hpp>

using stan::math::var;
using std::numeric_limits;
using std::vector;

class AgradCdfLogBeta : public AgradCdfLogTest {
 public:
  void valid_values(vector<vector<double> >& parameters,
                    vector<double>& log_cdf) {
    vector<double> param(3);

    param[0] = 0.5;  // y
    param[1] = 2.0;  // alpha (Success Scale)
    param[2] = 5.0;  // beta  (Failure Scale)
    parameters.push_back(param);
    log_cdf.push_back(std::log(0.890625));  // expected Log_CDF
  }

  void invalid_values(vector<size_t>& index, vector<double>& value) {
    // y
    index.push_back(0U);
    value.push_back(-1.0);

    index.push_back(0U);
    value.push_back(2.0);

    // alpha
    index.push_back(1U);
    value.push_back(-1.0);

    index.push_back(1U);
    value.push_back(0.0);

    index.push_back(1U);
    value.push_back(numeric_limits<double>::infinity());

    // beta
    index.push_back(2U);
    value.push_back(-1.0);

    index.push_back(2U);
    value.push_back(0.0);

    index.push_back(2U);
    value.push_back(numeric_limits<double>::infinity());
  }

  bool has_lower_bound() { return true; }

  double lower_bound() { return 0.0; }

  bool has_upper_bound() { return true; }

  double upper_bound() { return 1.0; }

  template <typename T_y, typename T_scale_succ, typename T_scale_fail,
            typename T3, typename T4, typename T5>
  stan::return_type_t<T_y, T_scale_succ, T_scale_fail> cdf_log(
      const T_y& y, const T_scale_succ& alpha, const T_scale_fail& beta,
      const T3&, const T4&, const T5&) {
    return stan::math::beta_cdf_log(y, alpha, beta);
  }

  template <typename T_y, typename T_scale_succ, typename T_scale_fail,
            typename T3, typename T4, typename T5>
  stan::return_type_t<T_y, T_scale_succ, T_scale_fail> cdf_log_function(
      const T_y& y, const T_scale_succ& alpha, const T_scale_fail& beta,
      const T3&, const T4&, const T5&) {
    return stan::math::beta_cdf_log(y, alpha, beta);
  }
};
