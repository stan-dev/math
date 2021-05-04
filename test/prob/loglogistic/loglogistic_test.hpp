// Arguments: Doubles, Doubles, Doubles
#include <stan/math/prim.hpp>

using stan::math::var;
using std::numeric_limits;
using std::vector;

class AgradDistributionsLoglogistic : public AgradDistributionTest {
 public:
  void valid_values(vector<vector<double> >& parameters,
                    vector<double>& log_prob) {
    vector<double> param(3);

    param[0] = 1.0;  // y
    param[1] = 2.0;  // alpha
    param[2] = 3.0;  // beta
    parameters.push_back(param);
    log_prob.push_back(-1.21639532432449293253);  // expected log_prob

    param[0] = 1.4;  // y
    param[1] = 1.6;  // alpha
    param[2] = 2.4;  // beta
    parameters.push_back(param);
    log_prob.push_back(-0.87286484103664152556);  // expected log_prob
  }

  void invalid_values(vector<size_t>& index, vector<double>& value) {
    // y
    index.push_back(0U);
    value.push_back(0.0);

    index.push_back(0U);
    value.push_back(-1.0);

    index.push_back(0U);
    value.push_back(numeric_limits<double>::infinity());

    index.push_back(0U);
    value.push_back(-numeric_limits<double>::infinity());

    // alpha
    index.push_back(1U);
    value.push_back(0.0);

    index.push_back(1U);
    value.push_back(-1.0);

    index.push_back(1U);
    value.push_back(numeric_limits<double>::infinity());

    index.push_back(1U);
    value.push_back(-numeric_limits<double>::infinity());

    // beta
    index.push_back(2U);
    value.push_back(0.0);

    index.push_back(2U);
    value.push_back(-1.0);

    index.push_back(2U);
    value.push_back(-numeric_limits<double>::infinity());

    index.push_back(2U);
    value.push_back(numeric_limits<double>::infinity());
  }

  template <typename T_y, typename T_scale, typename T_shape, typename T3,
            typename T4, typename T5>
  stan::return_type_t<T_y, T_scale, T_shape> log_prob(const T_y& y,
                                                      const T_scale& alpha,
                                                      const T_shape& beta,
                                                      const T3&, const T4&,
                                                      const T5&) {
    return stan::math::loglogistic_lpdf(y, alpha, beta);
  }

  template <bool propto, typename T_y, typename T_scale, typename T_shape,
            typename T3, typename T4, typename T5>
  stan::return_type_t<T_y, T_scale, T_shape> log_prob(const T_y& y,
                                                      const T_scale& alpha,
                                                      const T_shape& beta,
                                                      const T3&, const T4&,
                                                      const T5&) {
    return stan::math::loglogistic_lpdf<propto>(y, alpha, beta);
  }

  template <typename T_y, typename T_scale, typename T_shape, typename T3,
            typename T4, typename T5>
  stan::return_type_t<T_y, T_scale, T_shape> log_prob_function(
      const T_y& y, const T_scale& alpha, const T_shape& beta, const T3&,
      const T4&, const T5&) {
    using stan::math::log1p;
    return log(beta) - log(alpha) + (beta - 1) * (log(y) - log(alpha))
           - 2.0 * log1p(pow((y * stan::math::inv(alpha)), beta));
  }
};
