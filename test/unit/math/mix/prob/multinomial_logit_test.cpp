#include <test/unit/math/test_ad.hpp>

TEST(mathMixScalFun, multinomialLogit) {
  std::vector<int> ns{0, 1, 2, 3};
  Eigen::VectorXd beta(4);
  beta << 0.1, 0.1, 0.5, 0.3;

  auto f = [&ns](const auto& b) {
    return stan::math::multinomial_logit_lpmf(ns, b);
  };

  stan::test::expect_ad(f, beta);
}
