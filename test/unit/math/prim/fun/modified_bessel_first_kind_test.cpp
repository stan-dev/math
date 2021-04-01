#include <stan/math/prim.hpp>
#include <test/unit/math/prim/fun/binary_scalar_tester.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(MathFunctions, modified_bessel_first_kind) {
  using stan::math::modified_bessel_first_kind;

  EXPECT_FLOAT_EQ(11.301921952136330496356270183217102497412616594,
                  modified_bessel_first_kind(0, 4.0));
  EXPECT_FLOAT_EQ(-3.953370217402609396478635740580581287584221595,
                  modified_bessel_first_kind(1, -3.0));
  EXPECT_FLOAT_EQ(-3.953370217402609396478635740580581287584221595,
                  modified_bessel_first_kind(-1, -3.0));

  // compare integer argument to double argument
  EXPECT_FLOAT_EQ(modified_bessel_first_kind(0, 5),
                  modified_bessel_first_kind(0, 5.0));
  EXPECT_FLOAT_EQ(modified_bessel_first_kind(3, 1),
                  modified_bessel_first_kind(3, 1.0));
}

TEST(MathFunctions, modified_bessel_first_kind_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_THROW(stan::math::modified_bessel_first_kind(1, nan),
               std::domain_error);
}

TEST(MathFunctions, modified_bessel_first_kind_vec) {
  auto f = [](const auto& x1, const auto& x2) {
    using stan::math::modified_bessel_first_kind;
    return modified_bessel_first_kind(x1, x2);
  };

  std::vector<int> std_in1{1, 3, 1};
  Eigen::VectorXd in2(3);
  in2 << -1.3, 0.7, 2.8;
  stan::test::binary_scalar_tester(f, std_in1, in2);

  std::vector<std::vector<int>> std_std_in1{std_in1, std_in1, std_in1};
  Eigen::MatrixXd mat_in2 = in2.replicate(1, 3);
  stan::test::binary_scalar_tester(f, std_std_in1, mat_in2);
}
