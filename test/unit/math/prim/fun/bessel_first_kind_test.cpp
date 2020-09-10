#include <stan/math/prim.hpp>
#include <test/unit/math/prim/fun/binary_scalar_tester.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(MathFunctions, bessel_first_kind) {
  using stan::math::bessel_first_kind;

  EXPECT_FLOAT_EQ(-0.39714980986384737228659076845169804197561868528938,
                  bessel_first_kind(0, 4.0));
  EXPECT_FLOAT_EQ(-0.33905895852593645892551459720647889697308041819800,
                  bessel_first_kind(1, -3.0));
  EXPECT_FLOAT_EQ(0.33905895852593645892551459720647889697308041819800,
                  bessel_first_kind(-1, -3.0));
}

TEST(MathFunctions, bessel_first_kind_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_THROW(stan::math::bessel_first_kind(3, nan), std::domain_error);
}

TEST(MathFunctions, bessel_first_kind_vec) {
  auto f = [](const auto& x1, const auto& x2) {
    using stan::math::bessel_first_kind;
    return bessel_first_kind(x1, x2);
  };

  std::vector<int> std_in1{1, 3, 1};
  Eigen::VectorXd in2(3);
  in2 << -1.3, 0.7, 2.8;
  stan::test::binary_scalar_tester(f, std_in1, in2);

  Eigen::MatrixXd mat_in2 = in2.replicate(1, 3);
  std::vector<std::vector<int>> std_std_in1{std_in1, std_in1, std_in1};
  stan::test::binary_scalar_tester(f, std_std_in1, mat_in2);
}
