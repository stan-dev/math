#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>

using stan::math::check_nan_gpu;

// ---------- check_nan_gpu: matrix tests ----------
TEST(ErrorHandlingScalar,CheckFinite_Matrix) {
  const std::string function = "check_nan_gpu";
  Eigen::Matrix<double,Eigen::Dynamic,1> x;
  using stan::math::matrix_gpu;
  x.resize(3);
  x << -1, 0, 1;
  matrix_gpu xx(x);
  
  ASSERT_NO_THROW(check_nan_gpu(function, "xx", xx))
    << "check_nan_gpu should be true with finite xx";

  x.resize(3);
  x << -1, 0, std::numeric_limits<double>::infinity();
  matrix_gpu xx(x);
  EXPECT_THROW(check_nan_gpu(function, "xx", xx), std::domain_error)
    << "check_nan_gpu should throw exception on Inf";

  x.resize(3);
  x << -1, 0, -std::numeric_limits<double>::infinity();
  matrix_gpu xx(x);

  EXPECT_THROW(check_nan_gpu(function, "xx", xx), std::domain_error)
    << "check_nan_gpu should throw exception on -Inf";

  x.resize(3);
  x << -1, 0, std::numeric_limits<double>::quiet_NaN();
  matrix_gpu xx(x);
  EXPECT_THROW(check_nan_gpu(function, "xx", xx), std::domain_error)
    << "check_nan_gpu should throw exception on NaN";
}

TEST(ErrorHandlingScalar,CheckFinite_nan) {
  const std::string function = "check_nan_gpu";
  double nan = std::numeric_limits<double>::quiet_NaN();

  Eigen::Matrix<double,Eigen::Dynamic,1> x_mat(3);
  x_mat << nan, 0, 1;
  matrix_gpu xx_mat(x_mat);
  EXPECT_THROW(check_finite(function, "xx_mat", xx_mat),
               std::domain_error);

  x_mat << 1, nan, 1;
  matrix_gpu xx_mat(x_mat);
  EXPECT_THROW(check_finite(function, "xx_mat", xx_mat),
               std::domain_error);

  x_mat << 1, 0, nan;
  matrix_gpu xx_mat(x_mat);
  EXPECT_THROW(check_finite(function, "xx_mat", xx_mat),
               std::domain_error);
}

