// Arguments: Ints, Doubles, Doubles, Doubles
#include <stan/math/prim/prob/beta_neg_binomial_lcdf.hpp>
#include <stan/math/prim/fun/lbeta.hpp>
#include <stan/math/prim/fun/lgamma.hpp>

using stan::math::var;
using std::numeric_limits;
using std::vector;

class AgradCdfLogBetaNegBinomial : public AgradCdfLogTest {
 public:
  void valid_values(vector<vector<double>>& parameters,
                    vector<double>& cdf_log) {
    vector<double> param(4);

    param[0] = 0;    // n
    param[1] = 1.0;  // r
    param[2] = 5.0;  // alpha
    param[3] = 1.0;  // beta
    parameters.push_back(param);
    cdf_log.push_back(std::log(0.833333333333333));  // expected cdf_log
  }

  void invalid_values(vector<size_t>& index, vector<double>& value) {
    // n

    // r
    index.push_back(1U);
    value.push_back(0.0);

    index.push_back(1U);
    value.push_back(-1.0);

    index.push_back(1U);
    value.push_back(std::numeric_limits<double>::infinity());

    // alpha
    index.push_back(2U);
    value.push_back(0.0);

    index.push_back(2U);
    value.push_back(-1.0);

    index.push_back(2U);
    value.push_back(std::numeric_limits<double>::infinity());

    // beta
    index.push_back(3U);
    value.push_back(0.0);

    index.push_back(3U);
    value.push_back(-1.0);

    index.push_back(3U);
    value.push_back(std::numeric_limits<double>::infinity());
  }

  // BOUND INCLUDED IN ORDER FOR TEST TO PASS WITH CURRENT FRAMEWORK
  bool has_lower_bound() { return false; }

  bool has_upper_bound() { return false; }

  template <typename T_n, typename T_r, typename T_size1, typename T_size2,
            typename T4, typename T5>
  stan::return_type_t<T_r, T_size1, T_size2> cdf_log(const T_n& n, const T_r& r,
                                                     const T_size1& alpha,
                                                     const T_size2& beta,
                                                     const T4&, const T5&) {
    return stan::math::beta_neg_binomial_lcdf(n, r, alpha, beta);
  }

  template <typename T_n, typename T_r, typename T_size1, typename T_size2,
            typename T4, typename T5>
  stan::return_type_t<T_r, T_size1, T_size2> cdf_log_function(
      const T_n& n, const T_r& r, const T_size1& alpha, const T_size2& beta,
      const T4&, const T5&) {
    using stan::math::lbeta;
    using stan::math::lgamma;
    using stan::math::log_sum_exp;
    using std::vector;

    vector<stan::return_type_t<T_r, T_size1, T_size2>> lpmf_values;

    for (int i = 0; i <= n; i++) {
      auto lpmf = lbeta(i + r, alpha + beta) - lbeta(r, alpha)
                  + lgamma(i + beta) - lgamma(i + 1) - lgamma(beta);
      lpmf_values.push_back(lpmf);
    }

    return log_sum_exp(lpmf_values);
  }
};
