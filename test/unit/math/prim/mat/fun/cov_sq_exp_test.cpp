#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>

TEST(MathPrimMat, cov_sq_exp) {
  double sigma = 0.2;
  double l = 5;
  
  std::vector<double> x(3);
  x[0] = -2;
  x[1] = -1;
  x[2] = -0.5;


  Eigen::MatrixXd cov;
  EXPECT_NO_THROW(cov = stan::math::cov_sq_exp(x, sigma, l));

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      EXPECT_FLOAT_EQ(sigma * sigma * exp(pow(x[i] - x[j], 2) / (- 2.0 * l * l)),
                      cov(i, j))
        << "index: (" << i << ", " << j << ")";
}

TEST(MathPrimMat, error) {
  double sigma = 0.2;
  double l = 5;
  
  std::vector<double> x(3);
  x[0] = -2;
  x[1] = -1;
  x[2] = -0.5;

  double sigma_bad = -1;
  double l_bad = -1;
  std::vector<double> x_bad(1);
  x_bad[0] = std::numeric_limits<double>::quiet_NaN();

  EXPECT_THROW(stan::math::cov_sq_exp(x, sigma, l_bad), std::domain_error);
  EXPECT_THROW(stan::math::cov_sq_exp(x, sigma_bad, l), std::domain_error);
  EXPECT_THROW(stan::math::cov_sq_exp(x_bad, sigma, l), std::domain_error);
}

