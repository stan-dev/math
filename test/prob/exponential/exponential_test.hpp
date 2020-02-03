// Arguments: Doubles, Doubles
#include <stan/math/prim.hpp>

using stan::math::var;
using std::numeric_limits;
using std::vector;

class AgradDistributionsExponential : public AgradDistributionTest {
 public:
  void valid_values(vector<vector<double> >& parameters,
                    vector<double>& log_prob) {
    vector<double> param(2);

    param[0] = 2.0;  // y
    param[1] = 1.5;  // beta
    parameters.push_back(param);
    log_prob.push_back(-2.594534891891835393096);  // expected log_prob

    param[0] = 15.0;  // y
    param[1] = 3.9;   // beta
    parameters.push_back(param);
    log_prob.push_back(-57.13902344686439249699);  // expected log_prob

    param[0] = 1e-08;
    param[1] = 3.9;
    parameters.push_back(param);
    log_prob.push_back(1.3609765141356007);
  }

  void invalid_values(vector<size_t>& index, vector<double>& value) {
    // y
    index.push_back(0U);
    value.push_back(-10.0);

    index.push_back(0U);
    value.push_back(numeric_limits<double>::quiet_NaN());

    // beta
    index.push_back(1U);
    value.push_back(0.0);

    index.push_back(1U);
    value.push_back(-1.0);

    index.push_back(1U);
    value.push_back(numeric_limits<double>::infinity());

    index.push_back(1U);
    value.push_back(-numeric_limits<double>::infinity());
  }

  template <typename T_y, typename T_inv_scale, typename T2, typename T3,
            typename T4, typename T5>
  stan::return_type_t<T_y, T_inv_scale> log_prob(const T_y& y,
                                                 const T_inv_scale& beta,
                                                 const T2&, const T3&,
                                                 const T4&, const T5&) {
    return stan::math::exponential_log(y, beta);
  }

  template <bool propto, typename T_y, typename T_inv_scale, typename T2,
            typename T3, typename T4, typename T5>
  stan::return_type_t<T_y, T_inv_scale> log_prob(const T_y& y,
                                                 const T_inv_scale& beta,
                                                 const T2&, const T3&,
                                                 const T4&, const T5&) {
    return stan::math::exponential_log<propto>(y, beta);
  }

  template <typename T_y, typename T_inv_scale, typename T2, typename T3,
            typename T4, typename T5>
  stan::return_type_t<T_y, T_inv_scale> log_prob_function(
      const T_y& y, const T_inv_scale& beta, const T2&, const T3&, const T4&,
      const T5&) {
    return log(beta) - beta * y;
  }
};

TEST(ProbDistributionsExponential, Cumulative) {
  using stan::math::exponential_cdf;
  using std::numeric_limits;
  EXPECT_FLOAT_EQ(0.95021293, exponential_cdf(2.0, 1.5));
  EXPECT_FLOAT_EQ(1.0, exponential_cdf(15.0, 3.9));
  EXPECT_FLOAT_EQ(0.62280765, exponential_cdf(0.25, 3.9));

  // ??
  // EXPECT_FLOAT_EQ(0.0,
  //                 exponential_cdf(-numeric_limits<double>::infinity(),
  //                                 1.5));
  EXPECT_FLOAT_EQ(0.0, exponential_cdf(0.0, 1.5));
  EXPECT_FLOAT_EQ(1.0,
                  exponential_cdf(numeric_limits<double>::infinity(), 1.5));
}
