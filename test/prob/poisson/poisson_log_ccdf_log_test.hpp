// Arguments: Ints, Doubles
#include <stan/math/prim/prob/poisson_log_lccdf.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/lgamma.hpp>
#include <cmath>
#include <limits>

using stan::math::var;
using std::numeric_limits;
using std::vector;

class AgradCcdfLogPoissonLog : public AgradCcdfLogTest {
 public:
  void valid_values(vector<vector<double> >& parameters,
                    vector<double>& ccdf_log) {
    vector<double> param(2);

    param[0] = 17;              // n
    param[1] = std::log(13.0);  // alpha
    parameters.push_back(param);
    ccdf_log.push_back(
        std::log(1.0 - 0.8904649795242025600572));  // expected ccdf_log

    // ccdf_log_function below loses the digits of 1 - cdf where the ccdf is
    // small, so this point keeps the ccdf near 0.6
    param[0] = 40;              // n
    param[1] = std::log(42.0);  // alpha
    parameters.push_back(param);
    ccdf_log.push_back(-0.5413568944089692353016477);  // expected ccdf_log

    param[0] = 0.0;            // n
    param[1] = std::log(3.0);  // alpha
    parameters.push_back(param);
    ccdf_log.push_back(
        std::log(1.0 - 0.04978706836786394446248));  // expected ccdf_log
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
  stan::return_type_t<T_log_rate> ccdf_log(const T_n& n,
                                           const T_log_rate& alpha, const T2&,
                                           const T3&, const T4&, const T5&) {
    return stan::math::poisson_log_lccdf(n, alpha);
  }

  template <typename T_n, typename T_log_rate, typename T2, typename T3,
            typename T4, typename T5>
  stan::return_type_t<T_log_rate> ccdf_log_function(const T_n& n,
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
    return log(1.0 - cdf);
  }
};
