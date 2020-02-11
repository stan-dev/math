// Arguments: Ints, Ints, Ints
#include <stan/math/prim.hpp>

using std::vector;

class AgradDistributionsDiscreteRange : public AgradDistributionTest {
 public:
  void valid_values(vector<vector<double>>& parameters,
                    vector<double>& log_prob) {
    vector<double> param(3);

    param[0] = 3;   // y
    param[1] = 1;   // lower
    param[2] = 10;  // upper
    parameters.push_back(param);
    log_prob.push_back(-2.302585092994045901094);  // expected log_prob

    // case for lower == upper
    param[0] = 5;  // y
    param[1] = 5;  // lower
    param[2] = 5;  // upper
    parameters.push_back(param);
    log_prob.push_back(0);  // expected log_prob
  }

  void invalid_values(vector<size_t>& /*index*/, vector<double>& /*value*/) {
    // y

    // lower

    // upper
  }

  template <class T_y, class T_lower, class T_upper, class T3, typename T4,
            typename T5>
  stan::return_type_t<T_y, T_lower, T_upper> log_prob(const T_y& y,
                                                      const T_lower& lower,
                                                      const T_upper& upper,
                                                      const T3&, const T4&,
                                                      const T5&) {
    return stan::math::discrete_range_lpmf(y, lower, upper);
  }

  template <bool propto, class T_y, class T_lower, class T_upper, class T3,
            typename T4, typename T5>
  double log_prob(const T_y& y, const T_lower& lower, const T_upper& upper,
                  const T3&, const T4&, const T5&) {
    return stan::math::discrete_range_lpmf<propto>(y, lower, upper);
  }

  template <class T_y, class T_lower, class T_upper, class T3, typename T4,
            typename T5>
  double log_prob_function(const T_y& y, const T_lower& lower,
                           const T_upper& upper, const T3&, const T4&,
                           const T5&) {
    if (y < lower || y > upper) {
      return stan::math::LOG_ZERO;
    }

    return -log(upper - lower + 1);
  }
};
