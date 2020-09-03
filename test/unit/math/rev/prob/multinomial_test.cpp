#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>
#include <test/unit/math/rev/prob/expect_eq_diffs.hpp>
#include <gtest/gtest.h>
#include <string>
#include <vector>

template <typename T_prob>
void expect_propto_multinomial(std::vector<int>& ns1, T_prob theta1,
                               std::vector<int>& ns2, T_prob theta2,
                               std::string message) {
  expect_eq_diffs(stan::math::multinomial_log<false>(ns1, theta1),
                  stan::math::multinomial_log<false>(ns2, theta2),
                  stan::math::multinomial_log<true>(ns1, theta1),
                  stan::math::multinomial_log<true>(ns2, theta2), message);
}

TEST(AgradDistributionsMultinomial, Propto) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::var;
  std::vector<int> ns;
  ns.push_back(1);
  ns.push_back(2);
  ns.push_back(3);
  Matrix<var, Dynamic, 1> theta1(3, 1);
  theta1 << 0.3, 0.5, 0.2;
  Matrix<var, Dynamic, 1> theta2(3, 1);
  theta2 << 0.1, 0.2, 0.7;

  expect_propto_multinomial(ns, theta1, ns, theta2, "var: theta");
}

TEST(AgradDistributionsMultinomial, check_varis_on_stack) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::var;
  std::vector<int> ns;
  ns.push_back(1);
  ns.push_back(2);
  ns.push_back(3);
  Matrix<var, Dynamic, 1> theta(3, 1);
  theta << 0.3, 0.5, 0.2;

  test::check_varis_on_stack(stan::math::multinomial_log<false>(ns, theta));
  test::check_varis_on_stack(stan::math::multinomial_log<true>(ns, theta));
}
