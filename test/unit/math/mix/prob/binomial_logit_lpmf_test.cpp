#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>

TEST(mathMixScalFun, binomial_logit_lpmf) {
  auto f = [](const auto n, const auto N) {
    return [=](const auto& alpha) {
      return stan::math::binomial_logit_lpmf(n, N, alpha);
    };
  };

  Eigen::VectorXd alpha = Eigen::VectorXd::Random(3);
  std::vector<int> n_arr{1, 4, 5};
  std::vector<int> N_arr{10, 45, 25};

  stan::test::expect_ad(f(5, 25), 2.11);
  stan::test::expect_ad(f(5, 25), alpha);
  stan::test::expect_ad(f(n_arr, 25), alpha);
  stan::test::expect_ad(f(n_arr, N_arr), alpha);
  stan::test::expect_ad(f(n_arr, 10), 2.11);
  stan::test::expect_ad(f(n_arr, N_arr), 2.11);
  stan::test::expect_ad(f(5, N_arr), 2.11);
  stan::test::expect_ad(f(5, N_arr), alpha);
}
