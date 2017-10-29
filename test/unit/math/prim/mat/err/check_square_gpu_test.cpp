#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>

TEST(ErrorHandlingMatrix, check_square) {
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> y;
  
  y.resize(3,3);
  matrix_gpu yy(y);
  EXPECT_NO_THROW(stan::math::check_square("check_square",
                                           "yy", yy));

  y.resize(3, 2);
  matrix_gpu yy(y);
  EXPECT_THROW(stan::math::check_square("check_square", "yy", yy), 
               std::invalid_argument);
}

TEST(ErrorHandlingMatrix, check_square_nan) {
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> y;
  double nan = std::numeric_limits<double>::quiet_NaN();

  y.resize(3,3);
  y << nan, nan, nan,nan, nan, nan,nan, nan, nan;
  matrix_gpu yy(y);
  EXPECT_NO_THROW(stan::math::check_square("check_square",
                                           "yy", yy));

  y.resize(3, 2);
  y << nan, nan, nan,nan, nan, nan;
  matrix_gpu yy(y);
  EXPECT_THROW(stan::math::check_square("check_square", "yy", yy), 
               std::invalid_argument);
}

TEST(ErrorHandlingMatrix, check_square_0x0) {
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> y;
  
  y.resize(0,0);
  matrix_gpu yy(y);
  EXPECT_NO_THROW(stan::math::check_square("check_square",
                                           "yy", yy));
}

TEST(ErrorHandlingMatrix, check_square_0_size) {
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> y;
  
  y.resize(0,10);
  matrix_gpu yy(y);
  EXPECT_THROW(stan::math::check_square("check_square",
                                        "yy", yy),
               std::invalid_argument);

  y.resize(10,0);
  matrix_gpu yy(y);
  EXPECT_THROW(stan::math::check_square("check_square",
                                        "yy", yy),
               std::invalid_argument);

}
