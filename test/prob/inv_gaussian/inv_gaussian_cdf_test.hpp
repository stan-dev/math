// Arguments: Doubles, Doubles, Doubles
#include <stan/math/prim/prob/inv_gaussian_cdf.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/erfc.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/sqrt.hpp>

using stan::math::var;
using std::numeric_limits;
using std::vector;

class AgradCdfInvGaussian : public AgradCdfTest {
 public:
  void valid_values(vector<vector<double> >& parameters, vector<double>& cdf) {
    vector<double> param(3);

    param[0] = 1.2;  // y
    param[1] = 0.5;  // mu
    param[2] = 2.0;  // lambda
    parameters.push_back(param);
    cdf.push_back(0.9815922531042920910742);  // expected cdf

    param[0] = 0.3;  // y
    param[1] = 1.0;  // mu
    param[2] = 5.0;  // lambda
    parameters.push_back(param);
    cdf.push_back(0.003359190912064955560317);  // expected cdf
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

  bool has_lower_bound() { return true; }

  double lower_bound() { return 0.0; }

  bool has_upper_bound() { return false; }

  template <typename T_y, typename T_loc, typename T_shape, typename T3,
            typename T4, typename T5>
  stan::return_type_t<T_y, T_loc, T_shape> cdf(const T_y& y, const T_loc& mu,
                                               const T_shape& lambda, const T3&,
                                               const T4&, const T5&) {
    return stan::math::inv_gaussian_cdf(y, mu, lambda);
  }

  template <typename T_y, typename T_loc, typename T_shape, typename T3,
            typename T4, typename T5>
  stan::return_type_t<T_y, T_loc, T_shape> cdf_function(const T_y& y,
                                                        const T_loc& mu,
                                                        const T_shape& lambda,
                                                        const T3&, const T4&,
                                                        const T5&) {
    using stan::math::erfc;
    using stan::math::INV_SQRT_TWO;
    using std::exp;
    using std::sqrt;

    // Linear-space reference built from erfc, independent of the log-space
    // implementation under test.
    return 0.5 * erfc(-sqrt(lambda / y) * (y / mu - 1.0) * INV_SQRT_TWO)
           + exp(2.0 * lambda / mu) * 0.5
                 * erfc(sqrt(lambda / y) * (y / mu + 1.0) * INV_SQRT_TWO);
  }
};
