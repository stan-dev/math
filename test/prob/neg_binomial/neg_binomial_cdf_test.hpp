// Arguments: Ints, Doubles, Doubles
#include <stan/math/prim.hpp>
#include <boost/math/special_functions/binomial.hpp>

using stan::math::var;
using std::numeric_limits;
using std::vector;

class AgradCdfNegBinomial : public AgradCdfTest {
 public:
  void valid_values(vector<vector<double> >& parameters, vector<double>& cdf) {
    vector<double> param(3);

    param[0] = 15;  // Failures/Counts
    param[1] = 50;  // Successes/Shape
    param[2] = 3;   // logit(p)/Inverse Scale
    parameters.push_back(param);
    cdf.push_back(0.4240861277740262114122);  // expected cdf

    param[0] = 0;   // Failures/Counts
    param[1] = 15;  // Successes/Shape
    param[2] = 3;   // logit(p)/Inverse Scale
    parameters.push_back(param);
    cdf.push_back(0.013363461010158063716);  // expected cdf
  }

  void invalid_values(vector<size_t>& index, vector<double>& value) {
    // Successes/Shape
    index.push_back(1U);
    value.push_back(-1);

    // logit(p)/Inverse Scale
    index.push_back(2U);
    value.push_back(-1);
  }

  bool has_lower_bound() { return false; }

  bool has_upper_bound() { return false; }

  template <typename T_n, typename T_shape, typename T_inv_scale, typename T3,
            typename T4, typename T5>
  stan::return_type_t<T_shape, T_inv_scale> cdf(const T_n& n,
                                                const T_shape& alpha,
                                                const T_inv_scale& beta,
                                                const T3&, const T4&,
                                                const T5&) {
    return stan::math::neg_binomial_cdf(n, alpha, beta);
  }

  template <typename T_n, typename T_shape, typename T_inv_scale, typename T3,
            typename T4, typename T5>
  stan::return_type_t<T_shape, T_inv_scale> cdf_function(
      const T_n& n, const T_shape& alpha, const T_inv_scale& beta, const T3&,
      const T4&, const T5&) {
    using stan::math::binomial_coefficient_log;
    using std::exp;
    using std::log;

    stan::return_type_t<T_shape, T_inv_scale> cdf(0);

    for (int i = 0; i <= n; i++) {
      stan::return_type_t<T_shape, T_inv_scale> temp;
      temp = binomial_coefficient_log(i + alpha - 1, i);

      cdf += exp(temp + alpha * log(beta / (1 + beta))
                 + i * log(1 / (1 + beta)));
    }

    return cdf;
  }
};
