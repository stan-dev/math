// Arguments: Doubles, Doubles, Doubles
#include <stan/math/prim/prob/von_mises_lpdf.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/modified_bessel_first_kind.hpp>

using stan::math::var;
using std::numeric_limits;
using std::vector;

class AgradDistributionVonMises : public AgradDistributionTest {
 public:
  void valid_values(vector<vector<double> >& parameters,
                    vector<double>& log_prob) {
    vector<double> param(3);

    param[0] = boost::math::constants::third_pi<double>();  // y
    param[1] = boost::math::constants::sixth_pi<double>();  // mu
    param[2] = 0.5;                                         // kappa
    parameters.push_back(param);
    log_prob.push_back(-1.46641408370260739602);

    param[0] = -boost::math::constants::sixth_pi<double>();
    param[1] = -boost::math::constants::three_quarters_pi<double>();
    param[2] = 1.0;
    parameters.push_back(param);
    log_prob.push_back(-2.332610470019044868195);

    param[0] = boost::math::constants::pi<double>() / 4.;
    param[1] = -boost::math::constants::three_quarters_pi<double>();
    param[2] = 1.5;
    parameters.push_back(param);
    log_prob.push_back(-3.836664434124961609029);

    param[0] = -boost::math::constants::sixth_pi<double>();
    param[1] = boost::math::constants::sixth_pi<double>();
    param[2] = 4.0;
    parameters.push_back(param);
    log_prob.push_back(-2.262849861924804084623);

    param[0] = boost::math::constants::third_pi<double>();
    param[1] = boost::math::constants::sixth_pi<double>();
    param[2] = 1e-8;
    parameters.push_back(param);
    log_prob.push_back(-1.837877066409345339082);
  }

  void invalid_values(vector<size_t>& index, vector<double>& value) {
    // y
    index.push_back(0U);
    value.push_back(-numeric_limits<double>::infinity());

    index.push_back(0U);
    value.push_back(numeric_limits<double>::infinity());

    // mu
    index.push_back(1U);
    value.push_back(-numeric_limits<double>::infinity());

    index.push_back(1U);
    value.push_back(numeric_limits<double>::infinity());

    // kappa
    index.push_back(2U);
    value.push_back(-1.0);
  }

  template <typename T_y, typename T_loc, typename T_scale, typename T3,
            typename T4, typename T5>
  stan::return_type_t<T_y, T_loc, T_scale> log_prob(const T_y& y,
                                                    const T_loc& mu,
                                                    const T_scale& kappa,
                                                    const T3&, const T4&,
                                                    const T5&) {
    return stan::math::von_mises_lpdf(y, mu, kappa);
  }

  template <bool propto, typename T_y, typename T_loc, typename T_scale,
            typename T3, typename T4, typename T5>
  stan::return_type_t<T_y, T_loc, T_scale> log_prob(const T_y& y,
                                                    const T_loc& mu,
                                                    const T_scale& kappa,
                                                    const T3&, const T4&,
                                                    const T5&) {
    return stan::math::von_mises_lpdf<propto>(y, mu, kappa);
  }

  template <typename T_y, typename T_loc, typename T_scale, typename T3,
            typename T4, typename T5>
  stan::return_type_t<T_y, T_loc, T_scale> log_prob_function(
      const T_y& y, const T_loc& mu, const T_scale& kappa, const T3&, const T4&,
      const T5&) {
    using stan::math::modified_bessel_first_kind;
    using stan::math::pi;
    using std::log;

    return -stan::math::LOG_TWO_PI - log(modified_bessel_first_kind(0, kappa))
           + kappa * cos(mu - y);
  }
};
