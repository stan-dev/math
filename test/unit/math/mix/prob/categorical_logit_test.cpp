#include <test/unit/math/test_ad.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(ProbDistributionsCategoricalLogit, ValsTest) {
  auto f = [](const auto& x1) {
    return [=](const auto& x2) {
      return stan::math::categorical_logit_log(x1, x2);
    };
  };

  Eigen::VectorXd theta(3);
  theta << 0.2, 0.2, 0.6;

  std::vector<int> int_arr{1, 2, 1, 3};
  std::vector<Eigen::VectorXd> theta_arr(4);
  theta_arr[0] = Eigen::VectorXd::Random(5);
  theta_arr[1] = Eigen::VectorXd::Random(5);
  theta_arr[2] = Eigen::VectorXd::Random(5);
  theta_arr[3] = Eigen::VectorXd::Random(5);

  stan::test::expect_ad(f(1), theta);
  stan::test::expect_ad(f(int_arr), theta);
  stan::test::expect_ad(f(1), theta_arr);
  stan::test::expect_ad(f(int_arr), theta_arr);
}

TEST(ProbDistributionsCategoricalLogit, ErrsTest) {
  auto f = [](const auto& x1) {
    return [=](const auto& x2) {
      return stan::math::categorical_logit_log(x1, x2);
    };
  };

  std::vector<int> int_arr{1, 2, 1, 3};
  std::vector<int> int_arr_empty;
  std::vector<Eigen::VectorXd> theta_arr(2);
  theta_arr[0] = Eigen::VectorXd::Random(5);
  theta_arr[1] = Eigen::VectorXd::Random(5);

  stan::test::expect_ad(f(int_arr), theta_arr);
  stan::test::expect_ad(f(int_arr_empty), theta_arr);
}
