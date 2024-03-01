// Arguments: Doubles, Doubles, Doubles
#include <stan/math/prim/prob/pareto_lpdf.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/multiply_log.hpp>

using stan::math::var;
using std::numeric_limits;
using std::vector;

class AgradDistributionsPareto : public AgradDistributionTest {
 public:
  void valid_values(vector<vector<double> >& parameters,
                    vector<double>& log_prob) {
    vector<double> param(3);

    param[0] = 1.5;  // y
    param[1] = 0.5;  // y_min
    param[2] = 2.0;  // alpha
    parameters.push_back(param);
    log_prob.push_back(-1.909542504884438329782);  // expected log_prob

    param[0] = 19.5;  // y
    param[1] = 0.15;  // y_min
    param[2] = 5.0;   // alpha
    parameters.push_back(param);
    log_prob.push_back(-25.69864880541351226384);  // expected log_prob
  }

  void invalid_values(vector<size_t>& index, vector<double>& value) {
    // y

    // y_min
    index.push_back(1U);
    value.push_back(0.0);

    index.push_back(1U);
    value.push_back(-1.0);

    index.push_back(1U);
    value.push_back(numeric_limits<double>::infinity());

    index.push_back(1U);
    value.push_back(-numeric_limits<double>::infinity());

    // alpha
    index.push_back(2U);
    value.push_back(0.0);

    index.push_back(2U);
    value.push_back(-1.0);

    index.push_back(2U);
    value.push_back(numeric_limits<double>::infinity());

    index.push_back(2U);
    value.push_back(-numeric_limits<double>::infinity());
  }

  template <class T_y, class T_scale, class T_shape, typename T3, typename T4,
            typename T5>
  stan::return_type_t<T_y, T_scale, T_shape> log_prob(const T_y& y,
                                                      const T_scale& y_min,
                                                      const T_shape& alpha,
                                                      const T3&, const T4&,
                                                      const T5&) {
    return stan::math::pareto_lpdf(y, y_min, alpha);
  }

  template <bool propto, class T_y, class T_scale, class T_shape, typename T3,
            typename T4, typename T5>
  stan::return_type_t<T_y, T_scale, T_shape> log_prob(const T_y& y,
                                                      const T_scale& y_min,
                                                      const T_shape& alpha,
                                                      const T3&, const T4&,
                                                      const T5&) {
    return stan::math::pareto_lpdf<propto>(y, y_min, alpha);
  }

  template <class T_y, class T_scale, class T_shape, typename T3, typename T4,
            typename T5>
  stan::return_type_t<T_y, T_scale, T_shape> log_prob_function(
      const T_y& y, const T_scale& y_min, const T_shape& alpha, const T3&,
      const T4&, const T5&) {
    using stan::math::LOG_ZERO;
    using stan::math::multiply_log;

    if (y < y_min)
      return LOG_ZERO;
    return log(alpha) + multiply_log(alpha, y_min)
           - multiply_log(alpha + 1.0, y);
  }
};

TEST(ProbDistributionsParetoCDF, Values) {
  EXPECT_FLOAT_EQ(0.60434447, stan::math::pareto_cdf(3.45, 2.89, 5.235));
}
