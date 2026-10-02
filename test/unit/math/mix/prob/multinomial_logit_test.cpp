#include <test/unit/math/test_ad.hpp>
#include <limits>
#include <vector>

TEST_F(AgradRev, mathMixScalFun_multinomialLogit) {
  std::vector<int> ns{0, 1, 2, 3};
  Eigen::VectorXd beta(4);
  beta << 0.1, 0.1, 0.5, 0.3;

  auto f = [&ns](const auto& b) {
    return stan::math::multinomial_logit_lpmf(ns, b);
  };

  stan::test::expect_ad(f, beta);
}

TEST_F(AgradRev, mathMixScalFun_multinomialLogit_var_infinity_throws) {
  using stan::math::var;
  std::vector<int> ns{1, 2, 3};
  // autodiff args must be finite
  Eigen::Matrix<var, Eigen::Dynamic, 1> beta(3);
  beta << -std::numeric_limits<double>::infinity(), 0.1, 0.5;
  EXPECT_THROW(stan::math::multinomial_logit_lpmf(ns, beta), std::domain_error);
}
