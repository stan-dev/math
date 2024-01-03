#include <test/unit/math/test_ad.hpp>
#include <limits>

TEST(mathMixScalFun, bernoulliCDF) {
  // bind integer arg because can't autodiff through
  auto f = [](const auto& x1) {
    return [=](const auto& x2) { return stan::math::bernoulli_cdf(x1, x2); };
  };
  stan::test::expect_ad(f(0), 0.1);
  stan::test::expect_ad(f(0), std::numeric_limits<double>::quiet_NaN());
  stan::test::expect_ad(f(1), 0.5);
  stan::test::expect_ad(f(1), std::numeric_limits<double>::quiet_NaN());
  stan::test::expect_ad(f(1), 0.2);

  std::vector<int> std_in1{0, 1};
  Eigen::VectorXd in2(2);
  in2 << 0.5, 0.9;

  stan::test::expect_ad(f(std_in1), 0.2);
  stan::test::expect_ad(f(std_in1), std::numeric_limits<double>::quiet_NaN());
  stan::test::expect_ad(f(1), in2);
  stan::test::expect_ad(f(std_in1), in2);
  stan::test::expect_ad(f(std_in1), std::numeric_limits<double>::quiet_NaN());
}
