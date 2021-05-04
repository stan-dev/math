#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>
#include <test/unit/math/rev/prob/expect_eq_diffs.hpp>
#include <gtest/gtest.h>
#include <string>
#include <vector>

template <typename T_prob>
void expect_propto_multinomial_logit_lpmf(std::vector<int>& ns1, T_prob beta1,
                                          std::vector<int>& ns2, T_prob beta2,
                                          std::string message) {
  expect_eq_diffs(stan::math::multinomial_logit_lpmf<false>(ns1, beta1),
                  stan::math::multinomial_logit_lpmf<false>(ns2, beta2),
                  stan::math::multinomial_logit_lpmf<true>(ns1, beta1),
                  stan::math::multinomial_logit_lpmf<true>(ns2, beta2),
                  message);
}

TEST(AgradDistributionsMultinomialLogit, Propto) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::var;
  std::vector<int> ns;
  ns.push_back(1);
  ns.push_back(2);
  ns.push_back(3);
  Matrix<var, Dynamic, 1> beta1(3, 1);
  beta1 << log(0.3), log(0.5), log(0.2);
  Matrix<var, Dynamic, 1> beta2(3, 1);
  beta2 << log(0.1), log(0.2), log(0.7);

  expect_propto_multinomial_logit_lpmf(ns, beta1, ns, beta2, "var: beta");
}

TEST(AgradDistributionsMultinomialLogit, check_varis_on_stack) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::var;
  std::vector<int> ns;
  ns.push_back(1);
  ns.push_back(2);
  ns.push_back(3);
  Matrix<var, Dynamic, 1> beta(3, 1);
  beta << log(0.3), log(0.5), log(0.2);

  test::check_varis_on_stack(
      stan::math::multinomial_logit_lpmf<false>(ns, beta));
  test::check_varis_on_stack(
      stan::math::multinomial_logit_lpmf<true>(ns, beta));
}
