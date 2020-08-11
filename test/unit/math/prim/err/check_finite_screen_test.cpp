#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <string>
#include <vector>

TEST(ErrorHandlingArr, CheckFiniteScreen_Vector) {
  using stan::math::check_finite_screen;
  std::vector<double> x = {-1, 0, 1};
  EXPECT_TRUE(check_finite_screen(x));

  x = {-1, 0, std::numeric_limits<double>::infinity()};
  EXPECT_TRUE(check_finite_screen(x));

  x = {-1, 0, -std::numeric_limits<double>::infinity()};
  EXPECT_TRUE(check_finite_screen(x));

  x = {-1, 0, std::numeric_limits<double>::quiet_NaN()};
  EXPECT_TRUE(check_finite_screen(x));
}

TEST(ErrorHandlingArr, CheckFiniteScreen_std_vector_std_vector) {
  using stan::math::check_finite_screen;
  std::vector<double> x = {-1, 0, 1};
  std::vector<std::vector<double>> xx = {x};
  EXPECT_TRUE(check_finite_screen(xx));

  x = {-1, 0, std::numeric_limits<double>::infinity()};
  xx = {x};
  EXPECT_TRUE(check_finite_screen(xx));

  x = {-1, 0, -std::numeric_limits<double>::infinity()};
  xx = {x};
  EXPECT_TRUE(check_finite_screen(xx));

  x = {-1, 0, std::numeric_limits<double>::quiet_NaN()};
  xx = {x};
  EXPECT_TRUE(check_finite_screen(xx));
}

TEST(ErrorHandlingArr, CheckFiniteScreen_nan) {
  using stan::math::check_finite_screen;
  double nan = std::numeric_limits<double>::quiet_NaN();

  std::vector<double> x = {nan, 0, 1};
  EXPECT_TRUE(check_finite_screen(x));

  x = {1, nan, 1};
  EXPECT_TRUE(check_finite_screen(x));

  x = {1, 0, nan};
  EXPECT_TRUE(check_finite_screen(x));
}

TEST(ErrorHandlingMat, CheckFiniteScreen_Matrix) {
  using stan::math::check_finite_screen;
  Eigen::Matrix<double, Eigen::Dynamic, 1> x;

  x.resize(3);
  x << -1, 0, 1;
  EXPECT_FALSE(check_finite_screen(x));

  EXPECT_FALSE(check_finite_screen(x.array()));

  EXPECT_FALSE(check_finite_screen(x.transpose()));

  x.resize(3);
  x << -1, 0, std::numeric_limits<double>::infinity();
  EXPECT_TRUE(check_finite_screen(x));

  x.resize(3);
  x << -1, 0, -std::numeric_limits<double>::infinity();
  EXPECT_TRUE(check_finite_screen(x));

  x.resize(3);
  x << -1, 0, std::numeric_limits<double>::quiet_NaN();
  EXPECT_TRUE(check_finite_screen(x));
}

TEST(ErrorHandlingMat, CheckFiniteScreen_std_vector_Matrix) {
  using stan::math::check_finite_screen;
  Eigen::Matrix<double, Eigen::Dynamic, 1> x;

  x.resize(3);
  x << -1, 0, 1;

  std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1>> xv = {x};
  std::vector<Eigen::Array<double, Eigen::Dynamic, 1>> xva = {x.array()};
  std::vector<Eigen::Matrix<double, 1, Eigen::Dynamic>> xvt = {x.transpose()};

  EXPECT_TRUE(check_finite_screen(xv));

  EXPECT_TRUE(check_finite_screen(xva));

  EXPECT_TRUE(check_finite_screen(xvt));

  x.resize(3);
  x << -1, 0, std::numeric_limits<double>::infinity();
  xv = {x};
  EXPECT_TRUE(check_finite_screen(xv));

  x.resize(3);
  x << -1, 0, -std::numeric_limits<double>::infinity();
  xv = {x};
  EXPECT_TRUE(check_finite_screen(xv));

  x.resize(3);
  x << -1, 0, std::numeric_limits<double>::quiet_NaN();
  xv = {x};
  EXPECT_TRUE(check_finite_screen(xv));
}

TEST(ErrorHandlingMat, CheckFiniteScreen_nan) {
  using stan::math::check_finite_screen;
  double nan = std::numeric_limits<double>::quiet_NaN();

  Eigen::Matrix<double, Eigen::Dynamic, 1> x_mat(3);
  x_mat << nan, 0, 1;
  EXPECT_TRUE(check_finite_screen(x_mat));

  x_mat << 1, nan, 1;
  EXPECT_TRUE(check_finite_screen(x_mat));

  x_mat << 1, 0, nan;
  EXPECT_TRUE(check_finite_screen(x_mat));
}

TEST(ErrorHandlingScalar, CheckFiniteScreen) {
  using stan::math::check_finite_screen;
  double x = 0;

  EXPECT_TRUE(check_finite_screen(x));

  x = std::numeric_limits<double>::infinity();
  EXPECT_TRUE(check_finite_screen(x));

  x = -std::numeric_limits<double>::infinity();
  EXPECT_TRUE(check_finite_screen(x));

  x = std::numeric_limits<double>::quiet_NaN();
  EXPECT_TRUE(check_finite_screen(x));
}

TEST(ErrorHandlingScalar, CheckFiniteScreen_nan) {
  using stan::math::check_finite_screen;
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_TRUE(check_finite_screen(nan));
}

TEST(ErrorHandlingScalar, CheckFiniteScreenVectorization) {
  using stan::math::check_finite_screen;
  Eigen::MatrixXd m = Eigen::MatrixXd::Constant(3, 2, 0);
  EXPECT_TRUE(
      check_finite_screen(std::vector<Eigen::MatrixXd>{m, m, m}));
  Eigen::MatrixXd m2 = m;
  m2(1, 1) = stan::math::INFTY;
  EXPECT_TRUE(
      check_finite_screen(std::vector<Eigen::MatrixXd>{m, m2, m}));
}
