// Arguments: Doubles, Doubles, Doubles
#include <stan/math/prim/prob/beta_proportion_lccdf.hpp>

using stan::math::var;
using std::numeric_limits;
using std::vector;

class AgradCcdfLogBetaProportion : public AgradCcdfLogTest {
 public:
  void valid_values(vector<vector<double> >& parameters,
                    vector<double>& log_ccdf) {
    vector<double> param(3);

    param[0] = 0.5;  // y
    param[1] = 0.4;  // mu (location)
    param[2] = 5.0;  // kappa (precision)
    parameters.push_back(param);
    log_ccdf.push_back(std::log(1.0 - 0.6875));  // expected Log_CDF

    param[0] = 0.25;  // y
    param[1] = 0.75;  // mu (location)
    param[2] = 1.4;   // kappa (precision)
    parameters.push_back(param);
    log_ccdf.push_back(
        std::log(1.0 - 0.08724396598527127));  // expected Log_CDF
  }

  void invalid_values(vector<size_t>& index, vector<double>& value) {
    // y
    index.push_back(0U);
    value.push_back(-1.0);

    index.push_back(0U);
    value.push_back(2.0);

    // mu
    index.push_back(1U);
    value.push_back(-1.0);

    index.push_back(1U);
    value.push_back(0.0);

    index.push_back(1U);
    value.push_back(1.0);

    index.push_back(1U);
    value.push_back(numeric_limits<double>::infinity());

    // kappa
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

  template <typename T_y, typename T_loc, typename T_prec, typename T3,
            typename T4, typename T5>
  stan::return_type_t<T_y, T_loc, T_prec> ccdf_log(const T_y& y,
                                                   const T_loc& mu,
                                                   const T_prec& kappa,
                                                   const T3&, const T4&,
                                                   const T5&) {
    return stan::math::beta_proportion_lccdf(y, mu, kappa);
  }

  template <typename T_y, typename T_loc, typename T_prec, typename T3,
            typename T4, typename T5>
  stan::return_type_t<T_y, T_loc, T_prec> ccdf_log_function(
      const T_y& y, const T_loc& mu, const T_prec& kappa, const T3&, const T4&,
      const T5&) {
    return stan::math::beta_proportion_lccdf(y, mu, kappa);
  }
};
