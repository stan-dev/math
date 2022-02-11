#include <stan/math/prim.hpp>
#include <test/unit/math/prim/fun/binary_scalar_tester.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

TEST(MathFunctions, modified_bessel_second_kind) {
  using stan::math::modified_bessel_second_kind;

  EXPECT_FLOAT_EQ(0.011159676085853024269745195979833489225,
                  modified_bessel_second_kind(0, 4.0));
  EXPECT_THROW(modified_bessel_second_kind(1, -3.0), std::domain_error);
  EXPECT_THROW(modified_bessel_second_kind(-1, -3.0), std::domain_error);
}

TEST(MathFunctions, modified_bessel_second_kind_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_TRUE(std::isnan(stan::math::modified_bessel_second_kind(0, nan)));
}

TEST(MathFunctions, modified_bessel_second_kind_vec) {
  auto f = [](const auto& x1, const auto& x2) {
    using stan::math::modified_bessel_second_kind;
    return modified_bessel_second_kind(x1, x2);
  };

  std::vector<int> std_in1{1, 3, 1};
  Eigen::VectorXd in2(3);
  in2 << 1.3, 0.7, 2.8;
  stan::test::binary_scalar_tester(f, std_in1, in2);

  std::vector<std::vector<int>> std_std_in1{std_in1, std_in1, std_in1};
  Eigen::MatrixXd mat_in2 = in2.replicate(1, 3);
  stan::test::binary_scalar_tester(f, std_std_in1, mat_in2);
}
