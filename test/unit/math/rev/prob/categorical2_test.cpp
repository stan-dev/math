#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/prob/expect_eq_diffs.hpp>
#include <string>

template <typename T_prob>
inline void expect_propto_categorical_log(unsigned int n1, T_prob theta1,
                                          unsigned int n2, T_prob theta2,
                                          std::string message) {
  expect_eq_diffs(stan::math::categorical_log<false>(n1, theta1),
                  stan::math::categorical_log<false>(n2, theta2),
                  stan::math::categorical_log<true>(n1, theta1),
                  stan::math::categorical_log<true>(n2, theta2), message);
}

TEST(AgradDistributionsCategorical, Propto) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::var;
  unsigned int n;
  Matrix<var, Dynamic, 1> theta1(3, 1);
  theta1 << 0.3, 0.5, 0.2;
  Matrix<var, Dynamic, 1> theta2(3, 1);
  theta2 << 0.1, 0.2, 0.7;

  n = 1;
  expect_propto_categorical_log(n, theta1, n, theta2, "var: theta");

  n = 2;
  expect_propto_categorical_log(n, theta1, n, theta2, "var: theta");

  n = 3;
  expect_propto_categorical_log(n, theta1, n, theta2, "var: theta");
}
