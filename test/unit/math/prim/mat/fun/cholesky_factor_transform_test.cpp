#include <stan/math/prim/mat.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

using Eigen::Matrix;
using Eigen::Dynamic;

TEST(ProbTransform,choleskyFactor) {
  using Eigen::Matrix;
  using Eigen::Dynamic;
  using stan::math::cholesky_factor_constrain;
  using stan::math::cholesky_factor_free;
  
  Matrix<double,Dynamic,1> x(3);
  x << 1, 2, 3;
  
  Matrix<double,Dynamic,Dynamic> y
    = cholesky_factor_constrain(x,2,2);

  Matrix<double,Dynamic,1> x2
    = cholesky_factor_free(y);
  
  EXPECT_EQ(x2.size(), x.size());
  EXPECT_EQ(x2.rows(), x.rows());
  EXPECT_EQ(x2.cols(), x.cols());
  for (int i = 0; i < 3; ++i)
    EXPECT_FLOAT_EQ(x(i), x2(i));
}
TEST(ProbTransform,choleskyFactorLogJacobian) {
  using Eigen::Matrix;
  using Eigen::Dynamic;
  using stan::math::cholesky_factor_constrain;

  double lp;
  Matrix<double,Dynamic,1> x(3);

  x.resize(1);
  x << 2.3;
  lp = 1.9;
  cholesky_factor_constrain(x,1,1,lp);
  EXPECT_FLOAT_EQ(1.9 + 2.3, lp);
  
  x.resize(3);
  x << 
    1, 
    2, 3;
  lp = 7.2;
  cholesky_factor_constrain(x,2,2,lp);
  EXPECT_FLOAT_EQ(7.2 + 1 + 3, lp);

  x.resize(6);
  x << 
    1.001,
    2, 3.01,
    4, 5, 6.1;
  lp = 1.2;
  cholesky_factor_constrain(x,3,3,lp);
  EXPECT_FLOAT_EQ(1.2 + 1.001 + 3.01 + 6.1, lp);

  x.resize(9);
  lp = 1.2;
  x << 
    1.001,
    2, 3.01,
    4, 5, 6.1,
    7, 8, 9;
  cholesky_factor_constrain(x,4,3,lp);
  EXPECT_FLOAT_EQ(1.2 + 1.001 + 3.01 + 6.1, lp);

}
TEST(ProbTransform,choleskyFactorConstrainError) {
  using Eigen::Matrix;
  using Eigen::Dynamic;
  using stan::math::cholesky_factor_constrain;

  Matrix<double,Dynamic,1> x(3);
  x << 1, 2, 3;
  EXPECT_THROW(cholesky_factor_constrain(x,9,9), std::invalid_argument);
  double lp = 0;
  EXPECT_THROW(cholesky_factor_constrain(x,9,9,lp), std::invalid_argument);
}
TEST(ProbTransform,choleskyFactorFreeError) {
  using Eigen::Matrix;
  using Eigen::Dynamic;
  using stan::math::cholesky_factor_free;

  Matrix<double,Dynamic,Dynamic> y(1,1);
  y.resize(1,1);
  y << -2;
  EXPECT_THROW(cholesky_factor_free(y), std::domain_error);

  y.resize(2,2);
  y << 1, 2, 3, 4;
  EXPECT_THROW(cholesky_factor_free(y), std::domain_error);

  y.resize(2,3);
  y << 1, 0, 0,
    2, 3, 0;
  EXPECT_THROW(cholesky_factor_free(y), std::domain_error);
}

