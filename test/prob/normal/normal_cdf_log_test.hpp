// Arguments: Doubles, Doubles, Doubles
#include <stan/math/prim.hpp>

using stan::math::square;
using stan::math::var;
using std::numeric_limits;
using std::vector;

class AgradCdfLogNormal : public AgradCdfLogTest {
 public:
  void valid_values(vector<vector<double> >& parameters,
                    vector<double>& cdf_log) {
    vector<double> param(3);

    param[0] = 0;  // y
    param[1] = 0;  // mu
    param[2] = 1;  // sigma
    parameters.push_back(param);
    cdf_log.push_back(-0.6931471805599452862268);  // expected cdf_log

    param[0] = 1;  // y
    param[1] = 0;  // mu
    param[2] = 1;  // sigma
    parameters.push_back(param);
    cdf_log.push_back(-0.1727537790234499326392);  // expected cdf_log

    param[0] = -2;   // y
    param[1] = 0;    // mu
    param[2] = 1.5;  // sigma
    parameters.push_back(param);
    cdf_log.push_back(-2.3945773661586438052);  // expected cdf_log

    param[0] = -3.5;  // y
    param[1] = 1.9;   // mu
    param[2] = 7.2;   // sigma
    parameters.push_back(param);
    cdf_log.push_back(-1.484448229919656192521);  // expected cdf_log

    param[0] = 7.5;  // y
    param[1] = 0;    // mu
    param[2] = 1;    // sigma
    parameters.push_back(param);
    cdf_log.push_back(-3.1908916729109475e-14);  // expected cdf_log
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
    value.push_back(-numeric_limits<double>::infinity());
  }

  bool has_lower_bound() { return false; }

  bool has_upper_bound() { return false; }

  template <typename T_y, typename T_loc, typename T_scale, typename T3,
            typename T4, typename T5>
  stan::return_type_t<T_y, T_loc, T_scale> cdf_log(const T_y& y,
                                                   const T_loc& mu,
                                                   const T_scale& sigma,
                                                   const T3&, const T4&,
                                                   const T5&) {
    return stan::math::normal_cdf_log(y, mu, sigma);
  }

  template <typename T_y, typename T_loc, typename T_scale, typename T3,
            typename T4, typename T5>
  stan::return_type_t<T_y, T_loc, T_scale> cdf_log_function(
      const T_y& y, const T_loc& mu, const T_scale& sigma, const T3&, const T4&,
      const T5&) {
    using stan::math::LOG_HALF;
    using stan::math::SQRT_PI;
    using stan::math::SQRT_TWO;
    using std::exp;
    using std::log;
    using std::log1p;
    using std::pow;

    stan::return_type_t<T_y, T_loc, T_scale> cdf_log(0.0);

    // define constants used in numerical approximation of log(CDF) for x<-20
    stan::return_type_t<T_y, T_loc, T_scale> p[6]
        = {-0.000658749161529837803157, -0.0160837851487422766278,
           -0.125781726111229246204,    -0.360344899949804439429,
           -0.305326634961232344035,    -0.0163153871373020978498};
    stan::return_type_t<T_y, T_loc, T_scale> q[6]
        = {0.00233520497626869185443, 0.0605183413124413191178,
           0.527905102951428412248,   1.87295284992346047209,
           2.56852019228982242072,    1.0};
    stan::return_type_t<T_y, T_loc, T_scale> scaled_diff(0.0);
    scaled_diff += (y - mu) / (sigma * SQRT_TWO);
    stan::return_type_t<T_y, T_loc, T_scale> sigma2pi(0.0);
    sigma2pi += sigma * SQRT_TWO;

    // use erfc() instead of erf() in order to retain precision since for x>0
    // erfc()->0
    if (scaled_diff > 0.0) {
      // CDF(x) = 1/2 + 1/2erf(x) = 1 - 1/2erfc(x)
      cdf_log += log1p(-0.5 * erfc(scaled_diff));
    } else if (scaled_diff > -20.0) {
      // CDF(x) = 1/2 - 1/2erf(-x) = 1/2erfc(-x)
      cdf_log += log(erfc(-scaled_diff)) + LOG_HALF;
    } else {
      // check whether scaled_diff^10 term will overflow
      if (10.0 * log(-scaled_diff) < log(std::numeric_limits<double>::max())) {
        // entering territory where erfc(-x)~0
        // need to use direct numerical approximation of cdf_log instead
        // the following based on W. J. Cody, Math. Comp. 23(107):631-638 (1969)
        // CDF(x) = 1/2erfc(-x)
        stan::return_type_t<T_y, T_loc, T_scale> temp_p = 0.0;
        stan::return_type_t<T_y, T_loc, T_scale> temp_q = 0.0;
        stan::return_type_t<T_y, T_loc, T_scale> x2 = square(scaled_diff);
        stan::return_type_t<T_y, T_loc, T_scale> pow_x = 0.0;
        for (int j = 0; j < 6; j++) {
          pow_x = pow(scaled_diff, -2 * j);
          temp_p = temp_p - p[j] * pow_x;
          temp_q = temp_q - q[j] * pow_x;
        }
        cdf_log += log(0.5) + log(1.0 / SQRT_PI + (temp_p / temp_q) / x2)
                   - log(-scaled_diff) - x2;
      } else {
        cdf_log += stan::math::negative_infinity();
      }
    }
    return cdf_log;
  }
};
