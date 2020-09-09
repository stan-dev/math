#include <stan/math/prim.hpp>
#include <test/unit/math/prim/fun/binary_scalar_tester.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

TEST(MathFunctions, bessel_second_kind) {
  using stan::math::bessel_second_kind;

  EXPECT_FLOAT_EQ(-0.01694073932506499190363513444715321824049258989801,
                  bessel_second_kind(0, 4.0));
  EXPECT_FLOAT_EQ(0.3246744247917999784370128392879532396692751433723549,
                  bessel_second_kind(1, 3.0));
  EXPECT_THROW(bessel_second_kind(-1, -3.0), std::domain_error);
}

TEST(MathFunctions, bessel_second_kind_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_TRUE(std::isnan(stan::math::bessel_second_kind(1, nan)));
}

TEST(MathFunctions, bessel_second_kind_vec) {
  auto f = [](const auto& x1, const auto& x2) {
    using stan::math::bessel_second_kind;
    return bessel_second_kind(x1, x2);
  };

  std::vector<int> std_in1{1, 3, 1};
  Eigen::VectorXd in2(3);
  in2 << 1.3, 0.7, 2.8;
  stan::test::binary_scalar_tester(f, std_in1, in2);

  std::vector<std::vector<int>> std_std_in1{std_in1, std_in1, std_in1};
  Eigen::MatrixXd mat_in2 = in2.replicate(1, 3);
  stan::test::binary_scalar_tester(f, std_std_in1, mat_in2);
}
