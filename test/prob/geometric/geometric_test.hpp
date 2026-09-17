// Arguments: Ints, Doubles
#include <stan/math/prim/prob/geometric_lpmf.hpp>
#include <stan/math/prim/fun/log1m.hpp>

using stan::math::var;
using std::numeric_limits;
using std::vector;

class AgradDistributionsGeometric : public AgradDistributionTest {
 public:
  void valid_values(vector<vector<double>>& parameters,
                    vector<double>& log_prob) {
    vector<double> param(2);

    // lpmf(n=0|theta=0.5) = log(0.5) = -0.693147...
    param[0] = 0;    // n
    param[1] = 0.5;  // theta
    parameters.push_back(param);
    log_prob.push_back(-0.69314718055994528623);

    // lpmf(n=2|theta=0.5) = log(0.5) + 2*log(0.5) = -2.079441...
    param[0] = 2;    // n
    param[1] = 0.5;  // theta
    parameters.push_back(param);
    log_prob.push_back(-2.07944154167983574766);

    // lpmf(n=4|theta=0.3) = log(0.3) + 4*log(0.7) = -2.630672...
    param[0] = 4;    // n
    param[1] = 0.3;  // theta
    parameters.push_back(param);
    log_prob.push_back(-2.63067258008086568566);
  }

  void invalid_values(vector<size_t>& index, vector<double>& value) {
    // n
    index.push_back(0U);
    value.push_back(-1);

    // theta
    index.push_back(1U);
    value.push_back(-0.1);

    index.push_back(1U);
    value.push_back(1.1);
  }

  template <class T_n, class T_prob, typename T2, typename T3, typename T4,
            typename T5>
  stan::return_type_t<T_prob> log_prob(const T_n& n, const T_prob& theta,
                                       const T2&, const T3&, const T4&,
                                       const T5&) {
    return stan::math::geometric_lpmf(n, theta);
  }

  template <bool propto, class T_n, class T_prob, typename T2, typename T3,
            typename T4, typename T5>
  stan::return_type_t<T_prob> log_prob(const T_n& n, const T_prob& theta,
                                       const T2&, const T3&, const T4&,
                                       const T5&) {
    return stan::math::geometric_lpmf<propto>(n, theta);
  }

  template <class T_n, class T_prob, typename T2, typename T3, typename T4,
            typename T5>
  stan::return_type_t<T_prob> log_prob_function(const T_n& n,
                                                const T_prob& theta, const T2&,
                                                const T3&, const T4&,
                                                const T5&) {
    using stan::math::log1m;
    using std::log;
    return log(theta) + n * log1m(theta);
  }
};
