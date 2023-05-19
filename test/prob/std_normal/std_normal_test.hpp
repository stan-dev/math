// Arguments: Doubles
#include <stan/math/prim/prob/std_normal_lpdf.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/sqrt.hpp>

using stan::math::var;
using std::numeric_limits;
using std::vector;

class AgradDistributionStdNormal : public AgradDistributionTest {
 public:
  void valid_values(vector<vector<double> >& parameters,
                    vector<double>& log_prob) {
    vector<double> param(1);

    param[0] = 0;  // y
    parameters.push_back(param);
    log_prob.push_back(-0.918938533204672669541);  // expected log_prob

    param[0] = 1;  // y
    parameters.push_back(param);
    log_prob.push_back(-1.418938533204672669541);  // expected log_prob

    param[0] = -2;  // y
    parameters.push_back(param);
    log_prob.push_back(-2.918938533204672669541);  // expected log_prob
  }

  void invalid_values(vector<size_t>& index, vector<double>& value) {}

  template <typename T_y, typename T1, typename T2, typename T3, typename T4,
            typename T5>
  stan::return_type_t<T_y> log_prob(const T_y& y, const T1&, const T2&,
                                    const T3&, const T4&, const T5&) {
    return stan::math::std_normal_lpdf(y);
  }

  template <bool propto, typename T_y, typename T1, typename T2, typename T3,
            typename T4, typename T5>
  stan::return_type_t<T_y> log_prob(const T_y& y, const T1&, const T2&,
                                    const T3&, const T4&, const T5&) {
    return stan::math::std_normal_lpdf<propto>(y);
  }

  template <typename T_y, typename T1, typename T2, typename T3, typename T4,
            typename T5>
  stan::return_type_t<T_y, T1, T2, T3, T4, T5> log_prob_function(
      const T_y& y, const T1&, const T2&, const T3&, const T4&, const T5&) {
    using stan::math::pi;
    return -0.5 * y * y - log(sqrt(2.0 * pi()));
  }
};
