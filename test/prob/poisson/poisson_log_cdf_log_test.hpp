// Arguments: Ints, Doubles
#include <stan/math/prim/prob/poisson_log_lcdf.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/lgamma.hpp>
#include <cmath>
#include <limits>

using stan::math::var;
using std::numeric_limits;
using std::vector;

class AgradCdfLogPoissonLog : public AgradCdfLogTest {
 public:
  void valid_values(vector<vector<double> >& parameters,
                    vector<double>& cdf_log) {
    vector<double> param(2);

    param[0] = 17;              // n
    param[1] = std::log(13.0);  // alpha
    parameters.push_back(param);
    cdf_log.push_back(std::log(0.8904649795242025600572));  // expected cdf_log

    param[0] = 82;              // n
    param[1] = std::log(42.0);  // alpha
    parameters.push_back(param);
    cdf_log.push_back(std::log(0.9999999845303266798879));  // expected cdf_log

    param[0] = 0.0;            // n
    param[1] = std::log(3.0);  // alpha
    parameters.push_back(param);
    cdf_log.push_back(std::log(0.04978706836786394446248));  // expected cdf_log
  }

  void invalid_values(vector<size_t>& index, vector<double>& value) {
    // alpha
    index.push_back(1U);
    value.push_back(numeric_limits<double>::quiet_NaN());
  }

  bool has_lower_bound() { return false; }

  bool has_upper_bound() { return false; }

  template <typename T_n, typename T_log_rate, typename T2, typename T3,
            typename T4, typename T5>
  stan::return_type_t<T_log_rate> cdf_log(const T_n& n, const T_log_rate& alpha,
                                          const T2&, const T3&, const T4&,
                                          const T5&) {
    return stan::math::poisson_log_lcdf(n, alpha);
  }

  template <typename T_n, typename T_log_rate, typename T2, typename T3,
            typename T4, typename T5>
  stan::return_type_t<T_log_rate> cdf_log_function(const T_n& n,
                                                   const T_log_rate& alpha,
                                                   const T2&, const T3&,
                                                   const T4&, const T5&) {
    using stan::math::exp;
    using stan::math::lgamma;
    using stan::math::log;
    stan::return_type_t<T_log_rate> cdf(0);
    for (int i = 0; i <= n; i++) {
      cdf += exp(i * alpha - lgamma(i + 1));
    }
    cdf *= exp(-exp(alpha));
    return log(cdf);
  }
};
