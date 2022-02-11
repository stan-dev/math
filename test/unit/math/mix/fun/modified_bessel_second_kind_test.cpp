#include <test/unit/math/test_ad.hpp>
#include <limits>

TEST(mathMixScalFun, modifiedBesselSecondKind) {
  auto f = [](const int x1) {
    return [=](const auto& x2) {
      return stan::math::modified_bessel_second_kind(x1, x2);
    };
  };
  stan::test::expect_ad(f(1), -3.0);  // error

  stan::test::expect_ad(f(1), 4.0);
  stan::test::expect_ad(f(2), 2.0);

  stan::test::expect_ad(f(1), std::numeric_limits<double>::quiet_NaN());
}

TEST(mathMixScalFun, modifiedBesselSecondKind_vec) {
  auto f = [](const auto& x1, const auto& x2) {
    using stan::math::modified_bessel_second_kind;
    return modified_bessel_second_kind(x1, x2);
  };

  std::vector<int> std_in1{3, 1};
  Eigen::VectorXd in2(2);
  in2 << 0.5, 3.4;
  stan::test::expect_ad_vectorized_binary(f, std_in1, in2);

  std::vector<std::vector<int>> std_std_in1{std_in1, std_in1};
  Eigen::MatrixXd mat_in2 = in2.replicate(1, 2);
  stan::test::expect_ad_vectorized_binary(f, std_std_in1, mat_in2);
}
