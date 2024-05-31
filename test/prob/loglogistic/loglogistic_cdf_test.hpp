// Arguments: Doubles, Doubles, Doubles
#include <stan/math/prim/prob/loglogistic_cdf.hpp>
#include <stan/math/prim/fun/pow.hpp>

using stan::math::pow;
using stan::math::var;
using std::numeric_limits;
using std::vector;

class AgradCdfLogistic : public AgradCdfTest {
 public:
  void valid_values(vector<vector<double> >& parameters, vector<double>& cdf) {
    vector<double> param(3);

    param[0] = 2.0;  // y
    param[1] = 3.0;  // Scale
    param[2] = 2.0;  // Shape
    parameters.push_back(param);
    cdf.push_back(0.30769230769230765388);  // expected cdf

    param[0] = 20.0;  // y
    param[1] = 11.5;  // Scale
    param[2] = 15.0;  // Shape
    parameters.push_back(param);
    cdf.push_back(0.99975173823523311167);  // expected cdf

    param[0] = 0.0001;  // y
    param[1] = 1.2;     // Scale
    param[2] = 1.0;     // Shape
    parameters.push_back(param);
    cdf.push_back(0.00008332638946754438);  // expected cdf
  }

  void invalid_values(vector<size_t>& index, vector<double>& value) {
    // y
    index.push_back(0U);
    value.push_back(-1.0);

    index.push_back(0U);
    value.push_back(-numeric_limits<double>::infinity());

    // alpha
    index.push_back(1U);
    value.push_back(-1);

    index.push_back(1U);
    value.push_back(0);

    index.push_back(1U);
    value.push_back(numeric_limits<double>::infinity());

    // beta
    index.push_back(2U);
    value.push_back(-1.0);

    index.push_back(2U);
    value.push_back(0.0);

    index.push_back(2U);
    value.push_back(numeric_limits<double>::infinity());
  }

  bool has_lower_bound() { return true; }

  double lower_bound() { return 0.0; }

  bool has_upper_bound() { return false; }

  template <typename T_y, typename T_scale, typename T_shape, typename T3,
            typename T4, typename T5>
  stan::return_type_t<T_y, T_scale, T_shape> cdf(const T_y& y,
                                                 const T_scale& alpha,
                                                 const T_shape& beta, const T3&,
                                                 const T4&, const T5&) {
    return stan::math::loglogistic_cdf(y, alpha, beta);
  }

  template <typename T_y, typename T_scale, typename T_shape, typename T3,
            typename T4, typename T5>
  stan::return_type_t<T_y, T_scale, T_shape> cdf_function(const T_y& y,
                                                          const T_scale& alpha,
                                                          const T_shape& beta,
                                                          const T3&, const T4&,
                                                          const T5&) {
    return 1.0 / (1.0 + pow(alpha / y, beta));
  }
};
