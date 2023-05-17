// Arguments: Doubles, Doubles, Doubles
#include <stan/math/prim/prob/von_mises_cdf.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/floor.hpp>
using stan::math::var;
using std::numeric_limits;
using std::vector;

class AgradCdfVonMises : public AgradCdfTest {
 public:
  void valid_values(vector<vector<double> >& parameters, vector<double>& cdf) {
    vector<double> param(3);
    param[0] = -1.5707463267948965;  // y
    param[1] = 6.283185307179586;    // mu
    param[2] = 1.0;                  // kappa
    parameters.push_back(param);
    cdf.push_back(0.1097601896880533);  // expected CDF

    param[0] = -1.5707463267948965;  // y
    param[1] = 0;                    // mu
    param[2] = 0.01;                 // kappa
    parameters.push_back(param);
    cdf.push_back(0.2484164302237636);  // expected CDF

    param[0] = 0.0;                // y
    param[1] = 6.283185307179586;  // mu
    param[2] = 0.1;                // kappa
    parameters.push_back(param);
    cdf.push_back(0.5);  // expected CDF
  }

  void invalid_values(vector<size_t>& index, vector<double>& value) {
    // y
    index.push_back(0U);
    value.push_back(-stan::math::pi() - 1.0);

    index.push_back(0U);
    value.push_back(stan::math::pi() + 1.0);

    // mu
    index.push_back(1U);
    value.push_back(-numeric_limits<double>::infinity());

    index.push_back(1U);
    value.push_back(numeric_limits<double>::infinity());

    // kappa
    index.push_back(2U);
    value.push_back(-1.0);

    index.push_back(2U);
    value.push_back(-numeric_limits<double>::infinity());
  }

  bool has_lower_bound() { return true; }

  double lower_bound() { return -stan::math::pi(); }

  bool has_upper_bound() { return true; }

  double upper_bound() { return stan::math::pi(); }

  template <typename T_y, typename T_loc, typename T_scale, typename T3,
            typename T4, typename T5>
  stan::return_type_t<T_y, T_loc, T_scale> cdf(const T_y& y, const T_loc& mu,
                                               const T_scale& kappa, const T3&,
                                               const T4&, const T5&) {
    return stan::math::von_mises_cdf(y, mu, kappa);
  }

  template <typename T_x, typename T_mu, typename T_k, typename T3, typename T4,
            typename T5>
  stan::return_type_t<T_x, T_mu, T_k> cdf_function(const T_x& x, const T_mu& mu,
                                                   const T_k& k, const T3&,
                                                   const T4&, const T5&) {
    using return_t = stan::return_type_t<T_x, T_mu, T_k>;
    using stan::math::internal::von_mises_cdf_centered;
    const double pi = stan::math::pi();

    // shift x so that mean is 0
    return_t x2 = x - mu;

    // x2 is on an interval (2*n*pi, (2*n + 1)*pi), move it to (-pi, pi)
    x2 += pi;
    const auto x_floor = floor(x2 / stan::math::TWO_PI);
    const auto x_moded = x2 - x_floor * stan::math::TWO_PI;
    x2 = x_moded - pi;

    // mu is on an interval (2*n*pi, (2*n + 1)*pi), move it to (-pi, pi)
    T_mu mu2;
    mu2 = mu + pi;
    const auto mu_floor = floor(mu2 / stan::math::TWO_PI);
    const auto mu_moded = mu2 - mu_floor * stan::math::TWO_PI;
    mu2 = mu_moded - pi;

    // find cut
    return_t cut, bound_val;
    if (mu2 < 0) {
      cut = mu2 + pi;
      bound_val = -pi - mu2;
    }
    if (mu2 >= 0) {
      cut = mu2 - pi;
      bound_val = pi - mu2;
    }

    return_t f_bound_val = von_mises_cdf_centered(bound_val, k);
    if (x <= cut) {
      return von_mises_cdf_centered(x2, k) - f_bound_val;
    } else {
      return von_mises_cdf_centered(x2, k) + 1 - f_bound_val;
    }
  }
};
