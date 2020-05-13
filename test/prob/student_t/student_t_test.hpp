// Arguments: Doubles, Doubles, Doubles, Doubles
#include <stan/math/prim.hpp>

using stan::math::var;
using std::numeric_limits;
using std::vector;

class AgradDistributionsStudentT : public AgradDistributionTest {
 public:
  void valid_values(vector<vector<double> >& parameters,
                    vector<double>& log_prob) {
    vector<double> param(4);

    param[0] = 1.0;  // y
    param[1] = 1.0;  // nu
    param[2] = 0.0;  // mu
    param[3] = 1.0;  // sigma
    parameters.push_back(param);
    log_prob.push_back(-1.837877066409345339082);  // expected log_prob

    param[0] = -3.0;  // y
    param[1] = 2.0;   // nu
    param[2] = 0.0;   // mu
    param[3] = 1.0;   // sigma
    parameters.push_back(param);
    log_prob.push_back(-3.596842909197555560041);  // expected log_prob
  }

  void invalid_values(vector<size_t>& index, vector<double>& value) {
    // y

    // nu
    index.push_back(1U);
    value.push_back(0.0);

    index.push_back(1U);
    value.push_back(-1.0);

    index.push_back(1U);
    value.push_back(numeric_limits<double>::infinity());

    index.push_back(1U);
    value.push_back(-numeric_limits<double>::infinity());

    // mu
    index.push_back(2U);
    value.push_back(numeric_limits<double>::infinity());

    index.push_back(2U);
    value.push_back(-numeric_limits<double>::infinity());

    // sigma
    index.push_back(3U);
    value.push_back(0.0);

    index.push_back(3U);
    value.push_back(-1.0);

    index.push_back(3U);
    value.push_back(numeric_limits<double>::infinity());

    index.push_back(3U);
    value.push_back(-numeric_limits<double>::infinity());
  }

  template <class T_y, class T_dof, class T_loc, class T_scale, typename T4,
            typename T5>
  stan::return_type_t<T_y, T_dof, T_loc, T_scale> log_prob(
      const T_y& y, const T_dof& nu, const T_loc& mu, const T_scale& sigma,
      const T4&, const T5&) {
    return stan::math::student_t_log(y, nu, mu, sigma);
  }

  template <bool propto, class T_y, class T_dof, class T_loc, class T_scale,
            typename T4, typename T5>
  stan::return_type_t<T_y, T_dof, T_loc, T_scale> log_prob(
      const T_y& y, const T_dof& nu, const T_loc& mu, const T_scale& sigma,
      const T4&, const T5&) {
    return stan::math::student_t_log<propto>(y, nu, mu, sigma);
  }

  template <class T_y, class T_dof, class T_loc, class T_scale, typename T4,
            typename T5>
  stan::return_type_t<T_y, T_dof, T_loc, T_scale> log_prob_function(
      const T_y& y, const T_dof& nu, const T_loc& mu, const T_scale& sigma,
      const T4&, const T5&) {
    using boost::math::lgamma;
    using stan::math::log1p;
    using stan::math::LOG_SQRT_PI;
    using stan::math::square;
    using std::log;

    return lgamma((nu + 1.0) / 2.0) - lgamma(nu / 2.0) - LOG_SQRT_PI
           - 0.5 * log(nu) - log(sigma)
           - ((nu + 1.0) / 2.0) * log1p(square(((y - mu) / sigma)) / nu);
  }
};
