// Arguments: Doubles, Doubles, Doubles
#include <stan/math/prim/prob/inv_gaussian_lpdf.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/sqrt.hpp>
#include <stan/math/prim/fun/square.hpp>

using stan::math::var;
using std::numeric_limits;
using std::vector;

class AgradDistributionsInvGaussian : public AgradDistributionTest {
 public:
  void valid_values(vector<vector<double> >& parameters,
                    vector<double>& log_prob) {
    vector<double> param(3);

    param[0] = 1.2;  // y
    param[1] = 0.5;  // mu
    param[2] = 2.0;  // lambda
    parameters.push_back(param);
    log_prob.push_back(-2.479180611448965157415);  // expected log_prob

    param[0] = 0.3;  // y
    param[1] = 1.0;  // mu
    param[2] = 5.0;  // lambda
    parameters.push_back(param);
    log_prob.push_back(-2.391593703832052124008);  // expected log_prob
  }

  void invalid_values(vector<size_t>& index, vector<double>& value) {
    // y
    index.push_back(0U);
    value.push_back(-1.0);

    index.push_back(0U);
    value.push_back(-numeric_limits<double>::infinity());

    // mu
    index.push_back(1U);
    value.push_back(0.0);

    index.push_back(1U);
    value.push_back(-1.0);

    index.push_back(1U);
    value.push_back(numeric_limits<double>::infinity());

    index.push_back(1U);
    value.push_back(-numeric_limits<double>::infinity());

    // lambda
    index.push_back(2U);
    value.push_back(0.0);

    index.push_back(2U);
    value.push_back(-1.0);

    index.push_back(2U);
    value.push_back(numeric_limits<double>::infinity());

    index.push_back(2U);
    value.push_back(-numeric_limits<double>::infinity());
  }

  template <typename T_y, typename T_loc, typename T_shape, typename T3,
            typename T4, typename T5>
  stan::return_type_t<T_y, T_loc, T_shape> log_prob(const T_y& y,
                                                    const T_loc& mu,
                                                    const T_shape& lambda,
                                                    const T3&, const T4&,
                                                    const T5&) {
    return stan::math::inv_gaussian_lpdf(y, mu, lambda);
  }

  template <bool propto, typename T_y, typename T_loc, typename T_shape,
            typename T3, typename T4, typename T5>
  stan::return_type_t<T_y, T_loc, T_shape> log_prob(const T_y& y,
                                                    const T_loc& mu,
                                                    const T_shape& lambda,
                                                    const T3&, const T4&,
                                                    const T5&) {
    return stan::math::inv_gaussian_lpdf<propto>(y, mu, lambda);
  }

  template <typename T_y, typename T_loc, typename T_shape, typename T3,
            typename T4, typename T5>
  stan::return_type_t<T_y, T_loc, T_shape> log_prob_function(
      const T_y& y, const T_loc& mu, const T_shape& lambda, const T3&,
      const T4&, const T5&) {
    using stan::math::pi;
    using stan::math::square;
    using std::log;
    using std::sqrt;

    return log(sqrt(lambda / (2.0 * pi() * y * y * y)))
           - lambda * square(y - mu) / (2.0 * square(mu) * y);
  }
};
