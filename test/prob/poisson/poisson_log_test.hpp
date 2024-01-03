// Arguments: Ints, Doubles
#include <stan/math/prim/prob/poisson_log_lpmf.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/lgamma.hpp>

using stan::math::var;
using std::numeric_limits;
using std::vector;

class AgradDistributionsPoisson : public AgradDistributionTest {
 public:
  void valid_values(vector<vector<double> >& parameters,
                    vector<double>& log_prob) {
    using std::log;

    vector<double> param(2);

    param[0] = 17;         // n
    param[1] = log(13.0);  // alpha
    parameters.push_back(param);
    log_prob.push_back(-2.900934373290765311282);  // expected log_prob

    param[0] = 192;        // n
    param[1] = log(42.0);  // alpha
    parameters.push_back(param);
    log_prob.push_back(-145.3546649655311853166);  // expected log_prob

    param[0] = 0;         // n
    param[1] = log(3.0);  // alpha
    parameters.push_back(param);
    log_prob.push_back(-3.0);  // expected log_prob
  }

  void invalid_values(vector<size_t>& index, vector<double>& value) {
    // n
    index.push_back(0U);
    value.push_back(-1);

    // alpha
    // all OK
  }

  template <class T_n, class T_rate, typename T2, typename T3, typename T4,
            typename T5>
  stan::return_type_t<T_n, T_rate> log_prob(const T_n& n, const T_rate& alpha,
                                            const T2&, const T3&, const T4&,
                                            const T5&) {
    return stan::math::poisson_log_lpmf(n, alpha);
  }

  template <bool propto, class T_n, class T_rate, typename T2, typename T3,
            typename T4, typename T5>
  stan::return_type_t<T_n, T_rate> log_prob(const T_n& n, const T_rate& alpha,
                                            const T2&, const T3&, const T4&,
                                            const T5&) {
    return stan::math::poisson_log_lpmf<propto>(n, alpha);
  }

  template <class T_n, class T_rate, typename T2, typename T3, typename T4,
            typename T5>
  stan::return_type_t<T_n, T_rate> log_prob_function(const T_n& n,
                                                     const T_rate& alpha,
                                                     const T2&, const T3&,
                                                     const T4&, const T5&) {
    using boost::math::lgamma;
    using stan::math::LOG_ZERO;
    using stan::math::multiply_log;
    using std::exp;

    if (alpha == -std::numeric_limits<double>::infinity())
      return n == 0 ? 0 : LOG_ZERO;

    if (alpha == std::numeric_limits<double>::infinity())
      return LOG_ZERO;

    return -lgamma(n + 1.0) + n * alpha - exp(alpha);
  }
};
