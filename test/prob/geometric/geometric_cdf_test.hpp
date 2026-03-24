// Arguments: Ints, Doubles
#include <stan/math/prim/prob/geometric_cdf.hpp>
#include <stan/math/prim/fun/log1m.hpp>

using std::numeric_limits;
using std::vector;

class AgradCdfGeometric : public AgradCdfTest {
 public:
  void valid_values(vector<vector<double>>& parameters, vector<double>& cdf) {
    vector<double> param(2);

    // CDF(n=0|theta=0.5) = 1 - (0.5)^1 = 0.5
    param[0] = 0;    // n
    param[1] = 0.5;  // theta
    parameters.push_back(param);
    cdf.push_back(0.5);

    // CDF(n=2|theta=0.5) = 1 - (0.5)^3 = 0.875
    param[0] = 2;    // n
    param[1] = 0.5;  // theta
    parameters.push_back(param);
    cdf.push_back(0.875);

    // CDF(n=4|theta=0.3) = 1 - 0.7^5 = 0.83193
    param[0] = 4;    // n
    param[1] = 0.3;  // theta
    parameters.push_back(param);
    cdf.push_back(0.83193);

    // CDF(n=0|theta=1.0) = 1.0
    param[0] = 0;    // n
    param[1] = 1.0;  // theta
    parameters.push_back(param);
    cdf.push_back(1.0);
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
  stan::return_type_t<T_n, T_prob> cdf(const T_n& n, const T_prob& theta,
                                       const T2&, const T3&, const T4&,
                                       const T5&) {
    return stan::math::geometric_cdf(n, theta);
  }

  template <class T_n, class T_prob, typename T2, typename T3, typename T4,
            typename T5>
  stan::return_type_t<T_n, T_prob> cdf_function(const T_n& n,
                                                const T_prob& theta,
                                                const T2&, const T3&,
                                                const T4&, const T5&) {
    using stan::math::log1m;
    return 1.0 - stan::math::exp((n + 1.0) * log1m(theta));
  }
};
