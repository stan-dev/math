#include <stan/math/mix.hpp>
#include <test/unit/math/mix/prob/test_distribution_ad.hpp>

TEST(mathMixScalFun, normal_lpdf) {
  auto f = [](const double mu, const double sigma) {
    return [=](const auto& y) { return stan::math::normal_lpdf(y, mu, sigma); };
  };

  stan::test::expect_ad(f(0, 1), -2.3);
  stan::test::expect_ad(f(0, 1), 0.0);
  stan::test::expect_ad(f(0, 1), 1.7);
}

TEST(mathMixScalFun, vnormal_lpdf) {
  auto f = [](const auto& y, const auto& mu, const auto& sigma) {
    return stan::math::normal_lpdf<true, stan::math::ProbReturnType::Vector>(y, mu, sigma);
  };
  Eigen::VectorXd y = Eigen::VectorXd::Random(5);
  //y << 0, 0;
  Eigen::VectorXd mu = Eigen::VectorXd::Random(5);
  //mu << 0, 0;
  Eigen::VectorXd sigma = stan::math::abs(Eigen::VectorXd::Random(5));
  //sigma << 1, 1;
  stan::test::expect_ad_distribution(f, y, mu, sigma);
  //stan::test::expect_ad(f, y, mu, sigma);

}
