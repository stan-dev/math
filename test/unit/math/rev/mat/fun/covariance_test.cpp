#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/math/rev/mat/util.hpp>

TEST(AgradRevMatrix, cov_eq) {
  using stan::math::cov;
  using namespace Eigen;
    
  Matrix<double, Dynamic, Dynamic> x(3,3);
  x <<
    1, 5, 9,
    2, 4, 6,
    3, 6, 7;

  Matrix<double, Dynamic, Dynamic> cmat =
    cov(x);

  EXPECT_FLOAT_EQ(cmat(0,0), 1.0);
  EXPECT_FLOAT_EQ(cmat(0,2), -1.0);
  EXPECT_FLOAT_EQ(cmat(2,2), 7.0/3.0);
}


TEST(AgradRevMatrix, cov_colvec_eq) {
  using stan::math::cov;
  using namespace Eigen;  
  
  std::vector<Matrix<double, Dynamic,1> > mat_vec(3);
  for(int i = 0; i < 3; ++i)
    mat_vec[i] = Matrix<double, Dynamic, 1>(3);

  mat_vec[0] << 1, 2 ,3;
  mat_vec[1] << 5, 4, 6;
  mat_vec[2] << 9, 6, 7;
  
  Matrix<double, Dynamic, Dynamic> cmat =
    cov(mat_vec);

  EXPECT_FLOAT_EQ(cmat(0,0), 1.0);
  EXPECT_FLOAT_EQ(cmat(0,2), -1.0);
  EXPECT_FLOAT_EQ(cmat(2,2), 7.0/3.0);
}


TEST(AgradRevMatrix, cov_rowvec_eq) {
  using stan::math::cov;
  using namespace Eigen;
    
  std::vector<Matrix<double, 1, Dynamic> > mat_vec(3);
  for(int i = 0; i < 3; ++i)
    mat_vec[i] = Matrix<double,1,Dynamic>(3);

  mat_vec[0] << 1, 5 ,9;
  mat_vec[1] << 2, 4, 6;
  mat_vec[2] << 3, 6, 7;


  Matrix<double, Dynamic, Dynamic> cmat =
    cov(mat_vec);

  EXPECT_FLOAT_EQ(cmat(0,0), 1.0);
  EXPECT_FLOAT_EQ(cmat(0,2), -1.0);
  EXPECT_FLOAT_EQ(cmat(2,2), 7.0/3.0);
}
