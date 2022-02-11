#include <test/unit/math/test_ad.hpp>
#include <limits>

TEST(mathMixScalFun, lmgamma) {
  auto f = [](int x1) {
    return [=](const auto& x2) { return stan::math::hypot(x1, x2); };
  };
  stan::test::expect_ad(f(3), 3.2);
  stan::test::expect_ad(f(3), std::numeric_limits<double>::quiet_NaN());
}

TEST(mathMixScalFun, lmgamma_vec) {
  auto f = [](const auto& x1, const auto& x2) {
    using stan::math::lmgamma;
    return lmgamma(x1, x2);
  };

  std::vector<int> std_in1{3, 1};
  Eigen::VectorXd in2(2);
  in2 << 3.2, 3.4;
  stan::test::expect_ad_vectorized_binary(f, std_in1, in2);

  std::vector<std::vector<int>> std_std_in1{std_in1, std_in1};
  Eigen::MatrixXd mat_in2 = in2.replicate(1, 2);
  stan::test::expect_ad_vectorized_binary(f, std_std_in1, mat_in2);
}
