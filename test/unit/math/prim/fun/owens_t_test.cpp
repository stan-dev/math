#include <stan/math/prim.hpp>
#include <test/unit/math/prim/fun/binary_scalar_tester.hpp>
#include <boost/math/special_functions/owens_t.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

TEST(MathFunctions, owens_t) {
  double a = 1.0;
  double b = 2.0;
  EXPECT_FLOAT_EQ(stan::math::owens_t(a, b), boost::math::owens_t(a, b));
}

TEST(MathFunctions, owens_t_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_TRUE(std::isnan(stan::math::owens_t(1.0, nan)));
  EXPECT_TRUE(std::isnan(stan::math::owens_t(nan, 2.0)));
  EXPECT_TRUE(std::isnan(stan::math::owens_t(nan, nan)));
}

TEST(MathFunctions, owens_t_vec) {
  auto f = [](const auto& x1, const auto& x2) {
    using stan::math::owens_t;
    return owens_t(x1, x2);
  };

  Eigen::VectorXd in1(3);
  in1 << 1.8, 3.24, 1.8;
  Eigen::VectorXd in2(3);
  in2 << 7.3, 0.7, 2.8;
  stan::test::binary_scalar_tester(f, in1, in2);
}
