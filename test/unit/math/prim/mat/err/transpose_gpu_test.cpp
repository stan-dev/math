#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>

TEST(MathMatrix, transpose) {

  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> y;  
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> y_t;  
  y.resize(2,2);
  y << -1, 0, 1, 1;
  y_t << -1, 1, 0, 1;
  matrix_gpu yy(y);
  matrix_gpu yyy(2,2);
  using stan::math::transpose;
  EXPECT_NO_THROW(transpose(yyy,yy));
  EXPECT_NO_THROW(transpose(yy));
  stan::math::copy(yy, y);
  EXPECT_EQUAL(y, y_t);
  stan::math::copy(yyy, y);
  EXPECT_EQUAL(y, y_t);
  
}
