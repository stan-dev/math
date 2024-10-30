// Arguments: Ints, Doubles, Doubles, Doubles
#include <stan/math/prim/prob/beta_neg_binomial_lpmf.hpp>
#include <stan/math/prim/fun/lbeta.hpp>
#include <stan/math/prim/fun/lgamma.hpp>

using stan::math::var;
using std::numeric_limits;
using std::vector;

class AgradDistributionsBetaNegBinomial : public AgradDistributionTest {
 public:
  void valid_values(vector<vector<double> >& parameters,
                    vector<double>& log_prob) {
    vector<double> param(4);

    param[0] = 5;     // n
    param[1] = 20.0;  // r
    param[2] = 10.0;  // alpha
    param[3] = 25.0;  // beta
    parameters.push_back(param);
    log_prob.push_back(-10.3681267949788);  // expected log_prob

    param[0] = 10;   // n
    param[1] = 5.5;  // r
    param[2] = 2.5;  // alpha
    param[3] = 0.5;  // beta
    parameters.push_back(param);
    log_prob.push_back(-5.166741878823932);  // expected log_prob
  }

  void invalid_values(vector<size_t>& index, vector<double>& value) {
    // n
    index.push_back(0U);
    value.push_back(-1);

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

  template <class T_n, class T_r, class T_size1, class T_size2, typename T4,
            typename T5>
  stan::return_type_t<T_r, T_size1, T_size2> log_prob(const T_n& n,
                                                      const T_r& r,
                                                      const T_size1& alpha,
                                                      const T_size2& beta,
                                                      const T4&, const T5&) {
    return stan::math::beta_neg_binomial_lpmf(n, r, alpha, beta);
  }

  template <bool propto, class T_n, class T_r, class T_size1, class T_size2,
            typename T4, typename T5>
  stan::return_type_t<T_r, T_size1, T_size2> log_prob(const T_n& n,
                                                      const T_r& r,
                                                      const T_size1& alpha,
                                                      const T_size2& beta,
                                                      const T4&, const T5&) {
    return stan::math::beta_neg_binomial_lpmf<propto>(n, r, alpha, beta);
  }

  template <class T_n, class T_r, class T_size1, class T_size2, typename T4,
            typename T5>
  stan::return_type_t<T_r, T_size1, T_size2> log_prob_function(
      const T_n& n, const T_r& r, const T_size1& alpha, const T_size2& beta,
      const T4&, const T5&) {
    using stan::math::lbeta;
    using stan::math::lgamma;

    return lbeta(n + r, alpha + beta) - lbeta(r, alpha) + lgamma(n + beta)
           - lgamma(n + 1) - lgamma(beta);
  }
};
