// Arguments: Doubles, Doubles, Doubles, Doubles
#include <stan/math/prim/prob/generalized_normal_lpdf.hpp>
#include <stan/math/prim/fun/abs.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/pow.hpp>
#include <stan/math/prim/fun/tgamma.hpp>
#include <stan/math/prim/fun/constants.hpp>

using std::numeric_limits;
using std::vector;

class AgradDistributionGeneralizedNormal : public AgradDistributionTest {
 public:
  void valid_values(vector<vector<double> >& parameters,
                    vector<double>& log_prob) {
    vector<double> param(4);

    param[0] = 0;  // y
    param[1] = 0;  // mu
    param[2] = 1;  // alpha
    param[3] = 2;  // beta
    parameters.push_back(param);
    log_prob.push_back(
        -0.57236494292470008707171367567652935582);  // expected log_prob

    param[0] = 1;  // y
    param[1] = 0;  // mu
    param[2] = 1;  // alpha
    param[3] = 2;  // beta
    parameters.push_back(param);
    log_prob.push_back(
        -1.5723649429247000870717136756765293558);  // expected log_prob

    param[0] = -2;  // y
    param[1] = 0;   // mu
    param[2] = 1;   // alpha
    param[3] = 2;   // beta
    parameters.push_back(param);
    log_prob.push_back(
        -4.5723649429247000870717136756765293558);  // expected log_prob

    param[0] = -3.5;  // y
    param[1] = 1.9;   // mu
    param[2] = 7.2;   // alpha
    param[3] = 2;     // beta
    parameters.push_back(param);
    log_prob.push_back(
        -3.1089459689467097140959090592117462617);  // expected log_prob

    param[0] = 0;  // y
    param[1] = 0;  // mu
    param[2] = 1;  // alpha
    param[3] = 1;  // beta
    parameters.push_back(param);
    log_prob.push_back(
        -0.69314718055994530941723212145817656808);  // expected log_prob

    param[0] = 1;  // y
    param[1] = 0;  // mu
    param[2] = 1;  // alpha
    param[3] = 1;  // beta
    parameters.push_back(param);
    log_prob.push_back(
        -1.6931471805599453094172321214581765681);  // expected log_prob

    param[0] = -2;  // y
    param[1] = 0;   // mu
    param[2] = 1;   // alpha
    param[3] = 1;   // beta
    parameters.push_back(param);
    log_prob.push_back(
        -2.6931471805599453094172321214581765681);  // expected log_prob

    param[0] = -3.5;  // y
    param[1] = 1.9;   // mu
    param[2] = 7.2;   // alpha
    param[3] = 1;     // beta
    parameters.push_back(param);
    log_prob.push_back(
        -3.4172282065819549364414275049933934740);  // expected log_prob

    param[0] = 0;    // y
    param[1] = 0;    // mu
    param[2] = 1;    // alpha
    param[3] = 1.5;  // beta
    parameters.push_back(param);
    log_prob.push_back(
        -0.59083234759930449611508182336583846717);  // expected log_prob

    param[0] = 1;    // y
    param[1] = 0;    // mu
    param[2] = 1;    // alpha
    param[3] = 1.5;  // beta
    parameters.push_back(param);
    log_prob.push_back(
        -1.5908323475993044961150818233658384672);  // expected log_prob

    param[0] = -2;   // y
    param[1] = 0;    // mu
    param[2] = 1;    // alpha
    param[3] = 1.5;  // beta
    parameters.push_back(param);
    log_prob.push_back(
        -3.4192594723454945937184592717852346243);  // expected log_prob

    param[0] = -3.5;  // y
    param[1] = 1.9;   // mu
    param[2] = 7.2;   // alpha
    param[3] = 1.5;   // beta
    parameters.push_back(param);
    log_prob.push_back(
        -3.2144324264596431082120695849657575107);  // expected log_prob

    param[0] = 0.5;  // y
    param[1] = 0;    // mu
    param[2] = 1;    // alpha
    param[3] = 1;    // beta
    parameters.push_back(param);
    log_prob.push_back(
        -1.1931471805599453094172321214581765681);  // expected log_prob

    param[0] = 0.5;  // y
    param[1] = 0;    // mu
    param[2] = 1;    // alpha
    param[3] = 1.5;  // beta
    parameters.push_back(param);
    log_prob.push_back(
        -0.94438573819257864983001127257011830807);  // expected log_prob
  }

  void invalid_values(vector<size_t>& index, vector<double>& value) {
    // y

    // mu
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
    value.push_back(-numeric_limits<double>::infinity());

    // beta
    index.push_back(3U);
    value.push_back(0.0);

    index.push_back(3U);
    value.push_back(-1.0);

    index.push_back(3U);
    value.push_back(-numeric_limits<double>::infinity());
  }

  template <typename T_y, typename T_loc, typename T_scale, typename T_shape,
            typename T5, typename T6>
  stan::return_type_t<T_y, T_loc, T_scale, T_shape> log_prob(
      const T_y& y, const T_loc& mu, const T_scale& alpha, const T_shape& beta,
      const T5&, const T6&) {
    return stan::math::generalized_normal_lpdf(y, mu, alpha, beta);
  }

  template <bool propto, typename T_y, typename T_loc, typename T_scale,
            typename T_shape, typename T5, typename T6>
  stan::return_type_t<T_y, T_loc, T_scale, T_shape> log_prob(
      const T_y& y, const T_loc& mu, const T_scale& alpha, const T_shape& beta,
      const T5&, const T6&) {
    return stan::math::generalized_normal_lpdf<propto>(y, mu, alpha, beta);
  }

  template <typename T_y, typename T_loc, typename T_scale, typename T_shape,
            typename T5, typename T6>
  stan::return_type_t<T_y, T_loc, T_scale, T_shape> log_prob_function(
      const T_y& y, const T_loc& mu, const T_scale& alpha, const T_shape& beta,
      const T5&, const T6&) {
    using stan::math::abs;
    using stan::math::inv;
    using stan::math::lgamma;
    using stan::math::log;
    using stan::math::LOG_TWO;

    auto base = abs(y - mu) / alpha;
    bool at_zero = stan::math::value_of_rec(base) == 0;
    if (at_zero)
      base += 1;
    auto pow_term = pow(base, beta);
    if (at_zero)
      pow_term = 0;
    return -LOG_TWO - log(alpha) - lgamma(1.0 + inv(beta)) - pow_term;
  }
};

TEST(ProbDistributionsGeneralizedNormal, VectorWithYEqualsMu) {
  using Eigen::VectorXd;
  using stan::math::generalized_normal_lpdf;
  using stan::math::var;

  // y[1] == mu[1], other elements differ
  VectorXd y(3), mu(3);
  y << -1.0, 0.5, 2.0;
  mu << 0.0, 0.5, 0.0;
  double alpha = 1.0;
  double beta = 1.5;

  double lp = generalized_normal_lpdf(y, mu, alpha, beta);
  EXPECT_TRUE(std::isfinite(lp));

  // Same with var types to check autodiff
  std::vector<var> y_v(y.data(), y.data() + y.size());
  std::vector<var> mu_v(mu.data(), mu.data() + mu.size());
  var alpha_v = alpha;
  var beta_v = beta;
  var lp_v = generalized_normal_lpdf(y_v, mu_v, alpha_v, beta_v);
  lp_v.grad();
  for (size_t i = 0; i < y_v.size(); ++i) {
    EXPECT_TRUE(std::isfinite(y_v[i].adj()));
    EXPECT_TRUE(std::isfinite(mu_v[i].adj()));
  }
  EXPECT_TRUE(std::isfinite(alpha_v.adj()));
  EXPECT_TRUE(std::isfinite(beta_v.adj()));
}
