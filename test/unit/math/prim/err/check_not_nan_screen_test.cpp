#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <vector>
#include <string>
#include <stdexcept>

TEST(ErrorHandlingArr, CheckNotNanScreenVectorized) {
  using stan::math::check_not_nan_screen;
  int N = 5;
  std::vector<double> x(N);

  x.assign(N, 0);
  EXPECT_TRUE(check_not_nan_screen(x))
      << "check_not_nan_screen(vector) should be true with finite x: " << x[0];

  x.assign(N, std::numeric_limits<double>::infinity());
  EXPECT_TRUE(check_not_nan_screen(x))
      << "check_not_nan_screen(vector) should be true with x = Inf: " << x[0];

  x.assign(N, -std::numeric_limits<double>::infinity());
  EXPECT_TRUE(check_not_nan_screen(x))
      << "check_not_nan_screen(vector) should be true with x = -Inf: " << x[0];

  x.assign(N, std::numeric_limits<double>::quiet_NaN());
  EXPECT_TRUE(check_not_nan_screen(x))
      << "check_not_nan_screen(vector) should throw exception on NaN: " << x[0];
}

TEST(ErrorHandlingMatrix, checkNotNanScreenEigenRow) {
  using stan::math::check_not_nan_screen;
  stan::math::vector_d y;
  y.resize(3);
  y << 1, 2, 3;

  EXPECT_FALSE(stan::math::check_not_nan_screen(y));
  EXPECT_FALSE(stan::math::check_not_nan_screen(y));

  y(1) = std::numeric_limits<double>::quiet_NaN();
  EXPECT_TRUE(stan::math::check_not_nan_screen(y));
  EXPECT_TRUE(stan::math::check_not_nan_screen(y));
}

TEST(ErrorHandlingScalar, CheckNotNanScreen) {
  using stan::math::check_not_nan_screen;
  double x = 0;

  EXPECT_TRUE(check_not_nan_screen(x));

  x = std::numeric_limits<double>::infinity();
  EXPECT_TRUE(check_not_nan_screen(x));

  x = -std::numeric_limits<double>::infinity();
  EXPECT_TRUE(check_not_nan_screen(x));

  x = std::numeric_limits<double>::quiet_NaN();
  EXPECT_TRUE(check_not_nan_screen(x));
}

TEST(ErrorHandlingScalar, CheckNotNaNScreenVectorization) {
  using stan::math::check_not_nan_screen;
  Eigen::MatrixXd m = Eigen::MatrixXd::Constant(3, 2, 0);
  EXPECT_TRUE(
      check_not_nan_screen(std::vector<Eigen::MatrixXd>{m, m, m}));
  Eigen::MatrixXd m2 = m;
  m2(1, 1) = stan::math::NOT_A_NUMBER;
  EXPECT_TRUE(
      check_not_nan_screen(std::vector<Eigen::MatrixXd>{m, m2, m}));
}
