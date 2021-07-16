#include <test/unit/math/test_ad.hpp>

TEST(mathMixScalFun, ldexp) {
  auto f = [](const auto& x1) { return stan::math::ldexp(x1, 5); };

  stan::test::expect_ad(f, 3.1);
  stan::test::expect_ad(f, 0.0);
  stan::test::expect_ad(f, -1.5);
  stan::test::expect_ad(f, stan::math::INFTY);
  stan::test::expect_ad(f, stan::math::NOT_A_NUMBER);
}

TEST(mathMixScalFun, ldexp_vec) {
  auto f = [](const auto& x1, const auto& x2) {
    using stan::math::ldexp;
    return ldexp(x1, x2);
  };

  Eigen::VectorXd in1(2);
  in1 << 0.5, 3.4;
  std::vector<int> std_in2{3, 1};
  stan::test::expect_ad_vectorized_binary(f, in1, std_in2);

  Eigen::MatrixXd mat_in1 = in1.replicate(1, 2);
  std::vector<std::vector<int>> std_std_in2{std_in2, std_in2};
  stan::test::expect_ad_vectorized_binary(f, mat_in1, std_std_in2);
}
