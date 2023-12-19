#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>

TEST(mathMixScalFun, binomial_logit_glm_lpmf) {
  auto f = [](const auto n, const auto N) {
    return [=](const auto& x, const auto& alpha, const auto& beta) {
      return stan::math::binomial_logit_glm_lpmf(n, N, x, alpha, beta);
    };
  };

  std::vector<int> n_arr{1, 4};
  std::vector<int> N_arr{10, 45};
  Eigen::MatrixXd x = Eigen::MatrixXd::Random(2, 2);
  Eigen::RowVectorXd x_rowvec = x.row(0);
  Eigen::VectorXd alpha = Eigen::VectorXd::Random(2);
  Eigen::VectorXd beta = Eigen::VectorXd::Random(2);

  stan::test::expect_ad(f(n_arr[0], N_arr[0]), x, alpha, beta);
  stan::test::expect_ad(f(n_arr[0], N_arr), x, alpha, beta);
  stan::test::expect_ad(f(n_arr, N_arr[0]), x, alpha, beta);
  stan::test::expect_ad(f(n_arr, N_arr), x, alpha, beta);
  stan::test::expect_ad(f(n_arr[0], N_arr[0]), x, alpha[0], beta);
  stan::test::expect_ad(f(n_arr[0], N_arr), x, alpha[0], beta);
  stan::test::expect_ad(f(n_arr, N_arr[0]), x, alpha[0], beta);
  stan::test::expect_ad(f(n_arr, N_arr), x, alpha[0], beta);
  stan::test::expect_ad(f(n_arr[0], N_arr[0]), x_rowvec, alpha, beta);
  stan::test::expect_ad(f(n_arr[0], N_arr), x_rowvec, alpha, beta);
  stan::test::expect_ad(f(n_arr, N_arr[0]), x_rowvec, alpha, beta);
  stan::test::expect_ad(f(n_arr, N_arr), x_rowvec, alpha, beta);
  stan::test::expect_ad(f(n_arr[0], N_arr[0]), x_rowvec, alpha[0], beta);
  stan::test::expect_ad(f(n_arr[0], N_arr), x_rowvec, alpha[0], beta);
  stan::test::expect_ad(f(n_arr, N_arr[0]), x_rowvec, alpha[0], beta);
  stan::test::expect_ad(f(n_arr, N_arr), x_rowvec, alpha[0], beta);
}
