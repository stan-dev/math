// Arguments: Doubles, Doubles, Doubles
#include <stan/math/prim.hpp>

using stan::math::var;
using std::numeric_limits;
using std::vector;

class AgradDistributionsCauchy : public AgradDistributionTest {
 public:
  void valid_values(vector<vector<double> >& parameters,
                    vector<double>& log_prob) {
    vector<double> param(3);

    param[0] = 1.0;  // y
    param[1] = 0.0;  // mu
    param[2] = 1.0;  // sigma
    parameters.push_back(param);
    log_prob.push_back(-1.837877066409345339082);  // expected log_prob

    param[0] = -1.5;  // y
    param[1] = 0.0;   // mu
    param[2] = 1.0;   // sigma
    parameters.push_back(param);
    log_prob.push_back(-2.323384882191046330036);  // expected log_prob

    param[0] = -1.5;  // y
    param[1] = -1.0;  // mu
    param[2] = 1.0;   // sigma
    parameters.push_back(param);
    log_prob.push_back(-1.367873437163609873224);  // expected log_prob
  }

  void invalid_values(vector<size_t>& index, vector<double>& value) {
    // y

    // mu
    index.push_back(1U);
    value.push_back(numeric_limits<double>::infinity());

    index.push_back(1U);
    value.push_back(-numeric_limits<double>::infinity());

    // sigma
    index.push_back(2U);
    value.push_back(0.0);

    index.push_back(2U);
    value.push_back(-1.0);

    index.push_back(2U);
    value.push_back(numeric_limits<double>::infinity());

    index.push_back(2U);
    value.push_back(-numeric_limits<double>::infinity());
  }

  template <typename T_y, typename T_loc, typename T_scale, typename T3,
            typename T4, typename T5>
  typename stan::return_type<T_y, T_loc, T_scale>::type log_prob(
      const T_y& y, const T_loc& mu, const T_scale& sigma, const T3&, const T4&,
      const T5&) {
    return stan::math::cauchy_log(y, mu, sigma);
  }

  template <bool propto, typename T_y, typename T_loc, typename T_scale,
            typename T3, typename T4, typename T5>
  typename stan::return_type<T_y, T_loc, T_scale>::type log_prob(
      const T_y& y, const T_loc& mu, const T_scale& sigma, const T3&, const T4&,
      const T5&) {
    return stan::math::cauchy_log<propto>(y, mu, sigma);
  }

  template <typename T_y, typename T_loc, typename T_scale, typename T3,
            typename T4, typename T5>
  typename stan::return_type<T_y, T_loc, T_scale>::type log_prob_function(
      const T_y& y, const T_loc& mu, const T_scale& sigma, const T3&, const T4&,
      const T5&) {
    using stan::math::log1p;
    using stan::math::square;
    return -stan::math::LOG_PI - log(sigma) - log1p(square((y - mu) / sigma));
  }
};

TEST(ProbDistributionsCauchy, Cumulative) {
  using stan::math::cauchy_cdf;
  EXPECT_FLOAT_EQ(0.75, cauchy_cdf(1.0, 0.0, 1.0));
  EXPECT_FLOAT_EQ(0.187167, cauchy_cdf(-1.5, 0.0, 1.0));
  EXPECT_FLOAT_EQ(0.187167, cauchy_cdf(-2.5, -1.0, 1.0));
}
