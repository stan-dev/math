#include <test/unit/math/test_ad.hpp>
#include <cmath>
#include <limits>
#include <vector>

template <typename T>
void expect_arith_instantiate() {
  auto a1 = stan::math::pow(T(1.0), 1);
  auto b1 = stan::math::pow(T(1.0), 1.0);
  auto c1 = stan::math::pow(1, T(1.0));
  auto d1 = stan::math::pow(1.0, T(1.0));
  auto e1 = stan::math::pow(T(1.0), T(1.0));

  auto a2 = stan::math::pow(std::complex<T>(1.0), 1);
  auto b2 = stan::math::pow(std::complex<T>(1.0), 1.0);
  auto c2 = stan::math::pow(1, std::complex<T>(1.0));
  auto d2 = stan::math::pow(1.0, std::complex<T>(1.0));
  auto e2 = stan::math::pow(std::complex<T>(1.0), std::complex<T>(1.0));
  auto f2 = stan::math::pow(std::complex<double>(1.0), std::complex<T>(1.0));
  auto g2 = stan::math::pow(std::complex<T>(1.0), std::complex<double>(1.0));
}

// this one's been tricky to instantiate, so test all instances
TEST(mathMixScalFun, powInstantiations) {
  using stan::math::fvar;
  using stan::math::var;
  expect_arith_instantiate<double>();
  expect_arith_instantiate<var>();
  expect_arith_instantiate<fvar<double>>();
  expect_arith_instantiate<fvar<fvar<double>>>();
  expect_arith_instantiate<fvar<var>>();
  expect_arith_instantiate<fvar<fvar<var>>>();
}

TEST(mathMixScalFun, pow) {
  auto f
      = [](const auto& x1, const auto& x2) { return stan::math::pow(x1, x2); };

  stan::test::expect_ad(f, -0.4, 0.5);
  stan::test::expect_ad(f, 0.5, 0.5);
  stan::test::expect_ad(f, 0.5, 1.0);
  stan::test::expect_ad(f, 0.5, 1.2);
  stan::test::expect_ad(f, 0.5, 5.0);
  stan::test::expect_ad(f, 1.0, 2.0);
  stan::test::expect_ad(f, 3.0, 4.0);
  for (double y = -2; y <= 4; y += 0.5)
    stan::test::expect_ad(f, 4.0, y);

  double nan = std::numeric_limits<double>::quiet_NaN();
  stan::test::expect_ad(f, 1.0, nan);
  stan::test::expect_ad(f, nan, 1.0);
  stan::test::expect_ad(f, nan, nan);

  Eigen::VectorXd in1(3);
  in1 << 0.5, 3.4, 5.2;
  Eigen::VectorXd in2(3);
  in2 << 3.3, 0.9, 2.1;
  stan::test::expect_ad(f, in1, in2);
  stan::test::expect_ad(f, in1, 2.0);
  stan::test::expect_ad(f, 2.0, in1);

  stan::test::expect_ad_vectorized_binary(f, in1, in2);
}
