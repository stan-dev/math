#include <test/unit/math/test_ad.hpp>
#include <limits>

TEST(mathMixScalFun, modifiedBesselFirstKind) {
  // bind integer arg because can't autodiff through
  auto f = [](const int x1) {
    return [=](const auto& x2) {
      return stan::math::modified_bessel_first_kind(x1, x2);
    };
  };
  stan::test::expect_ad(f(-1), 3.0);
  stan::test::expect_ad(f(1), 4.0);
  stan::test::expect_ad(f(1), std::numeric_limits<double>::quiet_NaN());
  stan::test::expect_ad(f(2), 3.5);
  stan::test::expect_ad(f(2), std::numeric_limits<double>::quiet_NaN());

  // integer argument
  stan::test::expect_ad(f(0), 5);
  stan::test::expect_ad(f(1), 12);
  stan::test::expect_ad(f(8), 3);
}

TEST(mathMixScalFun, modifiedBesselFirstKind_vec) {
  auto f = [](const auto& x1, const auto& x2) {
    using stan::math::modified_bessel_first_kind;
    return modified_bessel_first_kind(x1, x2);
  };

  std::vector<int> std_in1{3, 1};
  Eigen::VectorXd in2(2);
  in2 << 0.5, 3.4;
  stan::test::expect_ad_vectorized_binary(f, std_in1, in2);

  std::vector<std::vector<int>> std_std_in1{std_in1, std_in1};
  Eigen::MatrixXd mat_in2 = in2.replicate(1, 2);
  stan::test::expect_ad_vectorized_binary(f, std_std_in1, mat_in2);
}
