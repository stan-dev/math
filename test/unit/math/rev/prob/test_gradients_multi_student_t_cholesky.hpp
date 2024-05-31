#ifndef TEST_UNIT_MATH_REV_PROB_TEST_GRADIENTS_MULTI_STUDENT_T_CHOLESKY
#define TEST_UNIT_MATH_REV_PROB_TEST_GRADIENTS_MULTI_STUDENT_T_CHOLESKY
#include <cmath>
#include <vector>
#include <iomanip>
#include <stdexcept>
#include <stan/math.hpp>

template <typename F, typename T_y, typename T_mu, typename T_sigma,
          typename T_nu>
std::vector<double> finite_diffs_multi_normal2(
    const F& fun, const std::vector<T_y>& vec_y,
    const std::vector<T_mu>& vec_mu, const std::vector<T_sigma>& vec_sigma,
    const T_nu& nu, double epsilon = 1e-6) {
  std::vector<double> diffs;
  diffs.reserve(vec_y.size() + vec_mu.size() + vec_sigma.size());

  std::vector<double> vec_y_plus = stan::math::value_of(vec_y);
  std::vector<double> vec_y_minus = vec_y_plus;
  std::vector<double> vec_mu_plus = stan::math::value_of(vec_mu);
  std::vector<double> vec_mu_minus = vec_mu_plus;
  std::vector<double> vec_sigma_plus = stan::math::value_of(vec_sigma);
  std::vector<double> vec_sigma_minus = vec_sigma_plus;
  double nu_plus = stan::math::value_of(nu);
  double nu_minus = nu_plus;

  if (!stan::is_constant_all<T_y>::value) {
    for (size_t i = 0; i < vec_y.size(); ++i) {
      double recover_vec_y_plus = vec_y_plus[i];
      double recover_vec_y_minus = vec_y_minus[i];
      vec_y_plus[i] += epsilon;
      vec_y_minus[i] -= epsilon;
      diffs.push_back(
          (fun(vec_y_plus, vec_mu_plus, vec_sigma_plus, nu_plus)
           - fun(vec_y_minus, vec_mu_minus, vec_sigma_minus, nu_minus))
          / (2 * epsilon));
      vec_y_plus[i] = recover_vec_y_plus;
      vec_y_minus[i] = recover_vec_y_minus;
    }
  }
  if (!stan::is_constant_all<T_mu>::value) {
    for (size_t i = 0; i < vec_mu.size(); ++i) {
      double recover_vec_mu_plus = vec_mu_plus[i];
      double recover_vec_mu_minus = vec_mu_minus[i];
      vec_mu_plus[i] += epsilon;
      vec_mu_minus[i] -= epsilon;
      diffs.push_back(
          (fun(vec_y_plus, vec_mu_plus, vec_sigma_plus, nu_plus)
           - fun(vec_y_minus, vec_mu_minus, vec_sigma_minus, nu_minus))
          / (2 * epsilon));
      vec_mu_plus[i] = recover_vec_mu_plus;
      vec_mu_minus[i] = recover_vec_mu_minus;
    }
  }
  if (!stan::is_constant_all<T_sigma>::value) {
    for (size_t i = 0; i < vec_sigma.size(); ++i) {
      double recover_vec_sigma_plus = vec_sigma_plus[i];
      double recover_vec_sigma_minus = vec_sigma_minus[i];
      vec_sigma_plus[i] += epsilon;
      vec_sigma_minus[i] -= epsilon;
      diffs.push_back(
          (fun(vec_y_plus, vec_mu_plus, vec_sigma_plus, nu_plus)
           - fun(vec_y_minus, vec_mu_minus, vec_sigma_minus, nu_minus))
          / (2 * epsilon));
      vec_sigma_plus[i] = recover_vec_sigma_plus;
      vec_sigma_minus[i] = recover_vec_sigma_minus;
    }
  }
  if (!stan::is_constant_all<T_nu>::value) {
    nu_plus += epsilon;
    nu_minus -= epsilon;
    diffs.push_back(
        (fun(vec_y_plus, vec_mu_plus, vec_sigma_plus, nu_plus)
         - fun(vec_y_minus, vec_mu_minus, vec_sigma_minus, nu_minus))
        / (2 * epsilon));
  }
  return diffs;
}

template <typename F, typename T_y, typename T_mu, typename T_sigma,
          typename T_nu>
std::vector<double> grad_multi_normal2(const F& fun,
                                       const std::vector<T_y>& vec_y,
                                       const std::vector<T_mu>& vec_mu,
                                       const std::vector<T_sigma>& vec_sigma,
                                       const T_nu& nu) {
  stan::math::var fx = fun(vec_y, vec_mu, vec_sigma, nu);
  std::vector<double> grad;
  std::vector<stan::math::var> vec_vars;
  if (!stan::is_constant_all<T_y>::value) {
    for (size_t i = 0; i < vec_y.size(); i++)
      vec_vars.push_back(vec_y[i]);
  }
  if (!stan::is_constant_all<T_mu>::value) {
    for (size_t i = 0; i < vec_mu.size(); i++)
      vec_vars.push_back(vec_mu[i]);
  }
  if (!stan::is_constant_all<T_sigma>::value) {
    for (size_t i = 0; i < vec_sigma.size(); i++)
      vec_vars.push_back(vec_sigma[i]);
  }
  if (!stan::is_constant_all<T_nu>::value)
    vec_vars.push_back(nu);

  fx.grad(vec_vars, grad);
  return grad;
}

template <typename F, typename T_y, typename T_mu, typename T_sigma,
          typename T_nu>
void test_grad_multi_student_t_cholesky(const F& fun,
                                        const std::vector<T_y>& vec_y,
                                        const std::vector<T_mu>& vec_mu,
                                        const std::vector<T_sigma>& vec_sigma,
                                        const T_nu& nu) {
  using std::fabs;
  std::vector<double> diffs_finite
      = finite_diffs_multi_normal2(fun, vec_y, vec_mu, vec_sigma, nu);
  std::vector<double> diffs_var
      = grad_multi_normal2(fun, vec_y, vec_mu, vec_sigma, nu);
  EXPECT_EQ(diffs_finite.size(), diffs_var.size());
  for (size_t i = 0; i < diffs_finite.size(); ++i) {
    double tolerance
        = 1e-6 * fmax(fabs(diffs_finite[i]), fabs(diffs_var[i])) + 1e-14;
    EXPECT_NEAR(diffs_finite[i], diffs_var[i], tolerance);
  }
}
#endif
