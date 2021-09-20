// Arguments: Doubles, Ints, Doubles, Doubles, Doubles, Doubles, Doubles
#include <stan/math/prim/prob/ddm_lcdf.hpp>

using std::vector;
using stan::math::INFTY;

class AgradCdfDdm : public AgradCdfTest {
public:
  void valid_values(vector<vector<double> >& parameters,
                    vector<double>& cdf) {
    vector<double> param(7);
    
    // each expected log_prob is calculated with the R package `fddm` as follows
    // fddm::pfddm(rt, response, a, v, t0, w, sv, log = TRUE)
    
    param[0] = 1.0;   // rt
    param[1] = 1;     // response
    param[2] = 1.0;   // a
    param[3] = -1.0;  // v
    param[4] = 0.0;   // t0
    param[5] = 0.5;   // w
    param[6] = 0.0;   // sv
    parameters.push_back(param);
    log_prob.push_back(-0.3189645693469165);  // expected log_prob
    
    param[0] = 1.0;   // rt
    param[1] = 2;     // response
    param[2] = 1.0;   // a
    param[3] = -1.0;  // v
    param[4] = 0.0;   // t0
    param[5] = 0.5;   // w
    param[6] = 0.0;   // sv
    parameters.push_back(param);
    log_prob.push_back(-1.318964569346917);  // expected log_prob
    
    param[0] = 1.0;   // rt
    param[1] = 1;     // response
    param[2] = 2.0;   // a
    param[3] = -1.0;  // v
    param[4] = 0.0;   // t0
    param[5] = 0.5;   // w
    param[6] = 0.0;   // sv
    parameters.push_back(param);
    log_prob.push_back(-0.4105356025317656);  // expected log_prob
    
    param[0] = 1.0;   // rt
    param[1] = 1;     // response
    param[2] = 1.0;   // a
    param[3] = 1.0;   // v
    param[4] = 0.0;   // t0
    param[5] = 0.5;   // w
    param[6] = 0.0;   // sv
    parameters.push_back(param);
    log_prob.push_back(-1.318964569346917);  // expected log_prob
    
    param[0] = 1.0;   // rt
    param[1] = 1;     // response
    param[2] = 1.0;   // a
    param[3] = -1.0;  // v
    param[4] = 0.5;   // t0
    param[5] = 0.5;   // w
    param[6] = 0.0;   // sv
    parameters.push_back(param);
    log_prob.push_back(-0.403297055469638);  // expected log_prob
    
    param[0] = 1.0;   // rt
    param[1] = 1;     // response
    param[2] = 1.0;   // a
    param[3] = -1.0;  // v
    param[4] = 0.0;   // t0
    param[5] = 0.2;   // w
    param[6] = 0.0;   // sv
    parameters.push_back(param);
    log_prob.push_back(-0.08206670817568271);  // expected log_prob
    
    param[0] = 1.0;   // rt
    param[1] = 1;     // response
    param[2] = 1.0;   // a
    param[3] = -1.0;  // v
    param[4] = 0.0;   // t0
    param[5] = 0.5;   // w
    param[6] = 1.0;   // sv
    parameters.push_back(param);
    log_prob.push_back(-0.3658772238413681);  // expected log_prob
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
  
  bool has_upper_bound() { return true; }
  
  double upper_bound() { return 0.0; }
  
  template <typename T_rt, typename T_response, typename T_a, typename T_v,
            typename T_t0, typename T_w, typename T_sv, typename T7>
  stan::return_type_t<T_rt, T_response, T_a, T_v, T_t0, T_w, T_sv> cdf(
      const T_rt& rt, const T_response& response, const T_a& a, const T_v& v,
      const T_t0& t0, const T_w& w, const T_sv& sv, const T7&) {
    return stan::math::ddm_lcdf(rt, response, a, v, t0, w, sv);
  }
  
  template <typename T_rt, typename T_response, typename T_a, typename T_v,
            typename T_t0, typename T_w, typename T_sv, typename T7>
  stan::return_type_t<T_rt, T_response, T_a, T_v, T_t0, T_w, T_sv, T7>
  cdf_function(const T_rt& rt, const T_response& response, const T_a& a,
               const T_v& v, const T_t0& t0, const T_w& w, const T_sv& sv,
               const T7&) {
    return stan::math::ddm_lcdf(rt, response, a, v, t0, w, sv);
  }
};
