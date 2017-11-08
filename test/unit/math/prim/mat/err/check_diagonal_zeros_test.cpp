#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>

TEST(ErrorHandlingMatrix, checkDiagZerosMatrix) {
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> y;
  
  y.resize(2,2);
  y << -1, 0, 1, 1;
  stan::math::matrix_gpu yy(y);
  EXPECT_NO_THROW(stan::math::check_diagonal_zeros("check_diagonal_zeros",
                                           "yy", yy));

  y << -1, 0, 0, 0;
  stan::math::matrix_gpu yyy(y);
  EXPECT_THROW(stan::math::check_diagonal_zeros("check_diagonal_zeros", "yyy", yyy), 
               std::domain_error);
}
