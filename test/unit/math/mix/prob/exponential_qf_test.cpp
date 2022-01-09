#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>

TEST(mathMixProb, exponential_qf) {
  auto f = [](const auto& p, const auto& beta) {
    return stan::math::exponential_qf(p, beta);
  };

  Eigen::VectorXd p(2);
  p << 0.2, 0.6;

  Eigen::VectorXd beta(2);
  beta << 16, 2;

  stan::test::expect_ad(f, p, beta);
  stan::test::expect_ad(f, 0.1, 8);
}
