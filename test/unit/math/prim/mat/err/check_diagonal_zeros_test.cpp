#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>

TEST(ErrorHandlingMatrix, checkDiagZerosMatrix) {
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> y;
  
  y.resize(2,2);
  y << -1, 0, 1, 1;
  matrix_gpu yy(y);
  EXPECT_NO_THROW(stan::math::check_diagonal_zeros("check_diagonal_zeros",
                                           "yy", yy));

  y << -1, 0, 0, 1;
  matrix_gpu yy(y);
  EXPECT_THROW(stan::math::check_diagonal_zeros("check_diagonal_zeros", "yy", yy), 
               std::domain_error);
}
