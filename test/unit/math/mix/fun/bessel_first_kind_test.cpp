#include <test/unit/math/test_ad.hpp>
#include <limits>

TEST(mathMixScalFun, besselFirstKind) {
  // bind integer arg because can't autodiff through
  auto f = [](const int x1) {
    return
        [=](const auto& x2) { return stan::math::bessel_first_kind(x1, x2); };
  };
  stan::test::expect_ad(f(0), 4.0);
  stan::test::expect_ad(f(0), std::numeric_limits<double>::quiet_NaN());
  stan::test::expect_ad(f(1), -3.0);
  stan::test::expect_ad(f(1), std::numeric_limits<double>::quiet_NaN());
  stan::test::expect_ad(f(2), 2.79);
}

TEST(mathMixScalFun, besselFirstKind_vec) {
  auto f = [](const auto& x1, const auto& x2) {
    using stan::math::bessel_first_kind;
    return bessel_first_kind(x1, x2);
  };

  Eigen::VectorXi in1(2);
  in1 << 3, 1;
  Eigen::VectorXd in2(2);
  in2 << 0.5, 3.4;
  stan::test::expect_ad_vectorized_binary(f, in1, in2);
}
