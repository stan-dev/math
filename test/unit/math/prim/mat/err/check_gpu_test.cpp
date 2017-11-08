#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>

using stan::math::check_nan_gpu;

// ---------- check_nan_gpu: matrix tests ----------
TEST(ErrorHandlingScalarGPU,Check_nan_Matrix) {
  const std::string function = "check_nan_gpu";
  Eigen::Matrix<double,Eigen::Dynamic,1> x;
  using stan::math::matrix_gpu;
  x.resize(3);
  x << -1, 0, 1;
  matrix_gpu xx(x);

  ASSERT_NO_THROW(check_nan_gpu(function, "xx", xx))
    << "check_nan_gpu should be true with finite xx";
}

TEST(ErrorHandlingScalarGPU,Check_nan_Matrix_quit_nan) {
  const std::string function = "check_nan_gpu";
  Eigen::Matrix<double,Eigen::Dynamic,1> x;
  using stan::math::matrix_gpu;

  x.resize(3);
  x << -1, 0, std::numeric_limits<double>::quiet_NaN();
  matrix_gpu xx(x);
  EXPECT_THROW(check_nan_gpu(function, "xx", xx), std::domain_error)
    << "check_nan_gpu should throw exception on NaN";
}

TEST(ErrorHandlingScalarGPU,Check_nan_positions) {
  const std::string function = "check_nan_gpu";
  double nan = std::numeric_limits<double>::quiet_NaN();
  using stan::math::matrix_gpu;
  Eigen::Matrix<double,Eigen::Dynamic,1> x_mat(3);
  x_mat << nan, 0, 1;
  matrix_gpu xx_mat1(x_mat);
  EXPECT_THROW(check_nan_gpu(function, "xx_mat1", xx_mat1),
               std::domain_error);

  x_mat << 1, nan, 1;
  matrix_gpu xx_mat2(x_mat);
  EXPECT_THROW(check_nan_gpu(function, "xx_mat2", xx_mat2),
               std::domain_error);

  x_mat << 1, 0, nan;
  matrix_gpu xx_mat3(x_mat);
  EXPECT_THROW(check_nan_gpu(function, "xx_mat3", xx_mat3),
               std::domain_error);
}

