#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/math/rev/mat/util.hpp>

TEST(AgradRevMatrix, cov_eq) {
  using stan::math::corr;
  using namespace Eigen;
    
  Matrix<double, Dynamic, Dynamic> x(3,3);
  x <<
    1, 5, 9,
    2, 4, 6,
    3, 6, 7;

  Matrix<double, Dynamic, Dynamic> cmat =
    corr(x);

  EXPECT_FLOAT_EQ(cmat(0,0), 1.0);
  EXPECT_FLOAT_EQ(cmat(0,2), -1.0 / sqrt(7.0 / 3.0 ));
  EXPECT_FLOAT_EQ(cmat(1,2), 0.5 / sqrt(7.0 / 3.0 ));
  
}

