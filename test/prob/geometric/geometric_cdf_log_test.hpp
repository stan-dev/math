// Arguments: Ints, Doubles
#include <stan/math/prim/prob/geometric_lcdf.hpp>
#include <stan/math/prim/fun/log1m.hpp>

using std::numeric_limits;
using std::vector;

class AgradCdfLogGeometric : public AgradCdfLogTest {
 public:
  void valid_values(vector<vector<double>>& parameters,
                    vector<double>& cdf_log) {
    vector<double> param(2);

    // lcdf(n=0|theta=0.5) = log(0.5) = -0.693147...
    param[0] = 0;    // n
    param[1] = 0.5;  // theta
    parameters.push_back(param);
    cdf_log.push_back(-0.69314718055994528623);

    // lcdf(n=2|theta=0.5) = log(0.875) = -0.133531...
    param[0] = 2;    // n
    param[1] = 0.5;  // theta
    parameters.push_back(param);
    cdf_log.push_back(-0.13353139262452262681);

    // lcdf(n=4|theta=0.3) = log(0.83193) = -0.184006...
    param[0] = 4;    // n
    param[1] = 0.3;  // theta
    parameters.push_back(param);
    cdf_log.push_back(-0.18400697631582843550);
  }

  void invalid_values(vector<size_t>& index, vector<double>& value) {
    // theta
    index.push_back(1U);
    value.push_back(-0.1);

    index.push_back(1U);
    value.push_back(1.1);
  }

  bool has_lower_bound() { return false; }

  bool has_upper_bound() { return false; }

  template <class T_n, class T_prob, typename T2, typename T3, typename T4,
            typename T5>
  stan::return_type_t<T_n, T_prob> cdf_log(const T_n& n, const T_prob& theta,
                                           const T2&, const T3&, const T4&,
                                           const T5&) {
    return stan::math::geometric_lcdf(n, theta);
  }

  template <class T_n, class T_prob, typename T2, typename T3, typename T4,
            typename T5>
  stan::return_type_t<T_n, T_prob> cdf_log_function(const T_n& n,
                                                    const T_prob& theta,
                                                    const T2&, const T3&,
                                                    const T4&, const T5&) {
    using stan::math::log1m;
    using stan::math::log1m_exp;
    return log1m_exp((n + 1.0) * log1m(theta));
  }
};
