// Arguments: Doubles, Ints, Doubles, Doubles, Doubles, Doubles, Doubles
#include <stan/math/prim/prob/ddm_lpdf.hpp>

using stan::math::INFTY;
using std::vector;

class AgradDistributionDdm : public AgradDistributionTest {
 public:
  void valid_values(vector<vector<double> >& parameters,
                    vector<double>& log_prob) {
    vector<double> param(7);

    // each expected log_prob is calculated with the R package `fddm` as follows
    // fddm::dfddm(rt, response, a, v, t0, w, sv, log = TRUE)

    param[0] = 1.0;   // rt
    param[1] = 1;     // response
    param[2] = 1.0;   // a
    param[3] = -1.0;  // v
    param[4] = 0.0;   // t0
    param[5] = 0.5;   // w
    param[6] = 0.0;   // sv
    parameters.push_back(param);
    log_prob.push_back(-3.790072391414288);  // expected log_prob

    param[0] = 1.0;   // rt
    param[1] = 2;     // response
    param[2] = 1.0;   // a
    param[3] = -1.0;  // v
    param[4] = 0.0;   // t0
    param[5] = 0.5;   // w
    param[6] = 0.0;   // sv
    parameters.push_back(param);
    log_prob.push_back(-4.790072391414288);  // expected log_prob

    param[0] = 1.0;   // rt
    param[1] = 1;     // response
    param[2] = 2.0;   // a
    param[3] = -1.0;  // v
    param[4] = 0.0;   // t0
    param[5] = 0.5;   // w
    param[6] = 0.0;   // sv
    parameters.push_back(param);
    log_prob.push_back(-0.9754202070046213);  // expected log_prob

    param[0] = 1.0;  // rt
    param[1] = 1;    // response
    param[2] = 1.0;  // a
    param[3] = 1.0;  // v
    param[4] = 0.0;  // t0
    param[5] = 0.5;  // w
    param[6] = 0.0;  // sv
    parameters.push_back(param);
    log_prob.push_back(-4.790072391414288);  // expected log_prob

    param[0] = 1.0;   // rt
    param[1] = 1;     // response
    param[2] = 1.0;   // a
    param[3] = -1.0;  // v
    param[4] = 0.5;   // t0
    param[5] = 0.5;   // w
    param[6] = 0.0;   // sv
    parameters.push_back(param);
    log_prob.push_back(-1.072671222447106);  // expected log_prob

    param[0] = 1.0;   // rt
    param[1] = 1;     // response
    param[2] = 1.0;   // a
    param[3] = -1.0;  // v
    param[4] = 0.0;   // t0
    param[5] = 0.2;   // w
    param[6] = 0.0;   // sv
    parameters.push_back(param);
    log_prob.push_back(-4.621465563321226);  // expected log_prob

    param[0] = 1.0;   // rt
    param[1] = 1;     // response
    param[2] = 1.0;   // a
    param[3] = -1.0;  // v
    param[4] = 0.0;   // t0
    param[5] = 0.5;   // w
    param[6] = 1.0;   // sv
    parameters.push_back(param);
    log_prob.push_back(-4.07414598169426);  // expected log_prob
  }

  void invalid_values(vector<size_t>& index, vector<double>& value) {
    // rt
    index.push_back(0U);
    value.push_back(0.0);

    index.push_back(0U);
    value.push_back(-1.0);

    index.push_back(0U);
    value.push_back(INFTY);

    index.push_back(0U);
    value.push_back(-INFTY);

    // response
    index.push_back(1U);
    value.push_back(0);

    index.push_back(1U);
    value.push_back(3);

    index.push_back(1U);
    value.push_back(-1);

    index.push_back(1U);
    value.push_back(INFTY);

    index.push_back(1U);
    value.push_back(-INFTY);

    // a
    index.push_back(2U);
    value.push_back(0.0);

    index.push_back(2U);
    value.push_back(-1.0);

    index.push_back(2U);
    value.push_back(INFTY);

    index.push_back(2U);
    value.push_back(-INFTY);

    // v
    index.push_back(3U);
    value.push_back(INFTY);

    index.push_back(3U);
    value.push_back(-INFTY);

    // t0
    index.push_back(4U);
    value.push_back(-1);

    index.push_back(4U);
    value.push_back(INFTY);

    index.push_back(4U);
    value.push_back(-INFTY);

    // w
    index.push_back(5U);
    value.push_back(-0.1);

    index.push_back(5U);
    value.push_back(0.0);

    index.push_back(5U);
    value.push_back(1.0);

    index.push_back(5U);
    value.push_back(1.1);

    index.push_back(5U);
    value.push_back(INFTY);

    index.push_back(5U);
    value.push_back(-INFTY);

    // sv
    index.push_back(6U);
    value.push_back(-1.0);

    index.push_back(6U);
    value.push_back(INFTY);

    index.push_back(6U);
    value.push_back(-INFTY);
  }

  template <typename T_rt, typename T_response, typename T_a, typename T_v,
            typename T_t0, typename T_w, typename T_sv, typename T7>
  stan::return_type_t<T_rt, T_response, T_a, T_v, T_t0, T_w, T_sv> log_prob(
      const T_rt& rt, const T_response& response, const T_a& a, const T_v& v,
      const T_t0& t0, const T_w& w, const T_sv& sv, const T7&) {
    return stan::math::ddm_lpdf(rt, response, a, v, t0, w, sv);
  }

  template <bool propto, typename T_rt, typename T_response, typename T_a,
            typename T_v, typename T_t0, typename T_w, typename T_sv,
            typename T7>
  stan::return_type_t<T_rt, T_response, T_a, T_v, T_t0, T_w, T_sv> log_prob(
      const T_rt& rt, const T_response& response, const T_a& a, const T_v& v,
      const T_t0& t0, const T_w& w, const T_sv& sv, const T7&) {
    return stan::math::ddm_lpdf<propto>(rt, response, a, v, t0, w, sv);
  }

  template <typename T_rt, typename T_response, typename T_a, typename T_v,
            typename T_t0, typename T_w, typename T_sv, typename T7>
  stan::return_type_t<T_rt, T_response, T_a, T_v, T_t0, T_w, T_sv, T7>
  log_prob_function(const T_rt& rt, const T_response& response, const T_a& a,
                    const T_v& v, const T_t0& t0, const T_w& w, const T_sv& sv,
                    const T7&) {
    return stan::math::ddm_lpdf<true>(rt, response, a, v, t0, w, sv);
  }
};
