// Arguments: Ints, Ints, Doubles, Doubles
#include <stan/math/prim/prob/beta_binomial_lpmf.hpp>
#include <stan/math/prim/fun/lbeta.hpp>
#include <boost/math/special_functions/binomial.hpp>

using stan::math::var;
using std::numeric_limits;
using std::vector;

class AgradDistributionsBetaBinomial : public AgradDistributionTest {
 public:
  void valid_values(vector<vector<double> >& parameters,
                    vector<double>& log_prob) {
    vector<double> param(4);

    param[0] = 5;     // n
    param[1] = 20;    // N
    param[2] = 10.0;  // alpha
    param[3] = 25.0;  // beta
    parameters.push_back(param);
    log_prob.push_back(-1.85400682167860331387425071679);  // expected log_prob

    param[0] = 25;    // n
    param[1] = 100;   // N
    param[2] = 30.0;  // alpha
    param[3] = 50.0;  // beta
    parameters.push_back(param);
    log_prob.push_back(
        -4.376695614336585231351487436303084);  // expected log_prob
  }

  void invalid_values(vector<size_t>& index, vector<double>& value) {
    // n

    // N
    index.push_back(1U);
    value.push_back(-1);

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

  template <class T_n, class T_N, class T_size1, class T_size2, typename T4,
            typename T5>
  stan::return_type_t<T_size1, T_size2> log_prob(const T_n& n, const T_N& N,
                                                 const T_size1& alpha,
                                                 const T_size2& beta, const T4&,
                                                 const T5&) {
    return stan::math::beta_binomial_lpmf(n, N, alpha, beta);
  }

  template <bool propto, class T_n, class T_N, class T_size1, class T_size2,
            typename T4, typename T5>
  stan::return_type_t<T_size1, T_size2> log_prob(const T_n& n, const T_N& N,
                                                 const T_size1& alpha,
                                                 const T_size2& beta, const T4&,
                                                 const T5&) {
    return stan::math::beta_binomial_lpmf<propto>(n, N, alpha, beta);
  }

  template <class T_n, class T_N, class T_size1, class T_size2, typename T4,
            typename T5>
  stan::return_type_t<T_size1, T_size2> log_prob_function(
      const T_n& n, const T_N& N, const T_size1& alpha, const T_size2& beta,
      const T4&, const T5&) {
    using stan::math::binomial_coefficient_log;
    using stan::math::lbeta;

    if (n < 0 || n > N)
      return 0.0;

    return binomial_coefficient_log(N, n) + lbeta(n + alpha, N - n + beta)
           - lbeta(alpha, beta);
  }
};

TEST(ProbDistributionsBetaBinomialCDF, Values) {
  EXPECT_FLOAT_EQ(0.8868204314,
                  stan::math::beta_binomial_cdf(49, 100, 1.349, 3.938));
}
