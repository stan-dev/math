// Arguments: Ints, Doubles
#include <stan/math/prim/prob/yule_simon_lccdf.hpp>
#include <stan/math/prim/fun/lgamma.hpp>

using std::numeric_limits;
using std::vector;

class AgradCcdfLogYuleSimon : public AgradCcdfLogTest {
 public:
  void valid_values(vector<vector<double>>& parameters,
                    vector<double>& cdf_log) {
    vector<double> param(2);

    param[0] = 5;     // n
    param[1] = 20.0;  // alpha
    parameters.push_back(param);
    cdf_log.push_back(std::log(1.0 - 0.9999811782420478));  // expected ccdf_log

    param[0] = 10;   // n
    param[1] = 5.5;  // alpha
    parameters.push_back(param);
    cdf_log.push_back(std::log(1.0 - 0.9997987132162779));  // expected ccdf_log
  }

  void invalid_values(vector<size_t>& index, vector<double>& value) {
    // n

    // alpha
    index.push_back(1U);
    value.push_back(0.0);

    index.push_back(1U);
    value.push_back(-1.0);

    index.push_back(1U);
    value.push_back(std::numeric_limits<double>::infinity());
  }

  // BOUND INCLUDED IN ORDER FOR TEST TO PASS WITH CURRENT FRAMEWORK
  bool has_lower_bound() { return false; }

  bool has_upper_bound() { return false; }

  template <class T_n, class T_alpha, typename T2, typename T3, typename T4,
            typename T5>
  stan::return_type_t<T_n, T_alpha> ccdf_log(const T_n& n, const T_alpha& alpha,
                                             const T2&, const T3&, const T4&,
                                             const T5&) {
    return stan::math::yule_simon_lccdf(n, alpha);
  }

  template <class T_n, class T_alpha, typename T2, typename T3, typename T4,
            typename T5>
  stan::return_type_t<T_n, T_alpha> ccdf_log_function(const T_n& n,
                                                      const T_alpha& alpha,
                                                      const T2&, const T3&,
                                                      const T4&, const T5&) {
    using stan::math::lgamma;

    auto log_ccdf
        = lgamma(alpha + 1.0) + lgamma(n + 1.0) - lgamma(n + alpha + 1.0);

    return log_ccdf;
  }
};
