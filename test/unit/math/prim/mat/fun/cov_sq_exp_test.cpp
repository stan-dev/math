#include <gtest/gtest.h>
#include <stan/math/prim/mat/fun/cov_sq_exp.hpp>



TEST(MathPrimMat, cov_sq_exp) {
  Eigen::MatrixXd x(3, 1);
  x(0) = -2;
  x(1) = -1;
  x(2) = -0.5;

  double sigma = 0.2;
  double l = 5;

  Eigen::MatrixXd cov;

  EXPECT_NO_THROW(cov = stan::math::cov_sq_exp(x, sigma, l));

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      EXPECT_FLOAT_EQ(sigma * sigma * exp(pow(x(i) - x(j), 2) / (- 2.0 * l * l)),
                      cov(i, j))
        << "index: (" << i << ", " << j << ")";
}

TEST(MathPrimMat, error) {
  
}
