#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <string>
#include <vector>


template<typename T_x, typename T_sigma>
std::string pull_msg(std::vector<T_x> x, T_sigma sigma) {
  std::string message;
  try {
    stan::math::cov_dot_prod(x, sigma);
  } catch (std::domain_error& e) {
    message = e.what();
  } catch (...) {
    message = "Threw the wrong exection";
  }
  return message;
}

template<typename T_x1, typename T_x2, typename T_sigma>
std::string pull_msg(std::vector<T_x1> x1, std::vector<T_x2> x2, T_sigma sigma) {
  std::string message;
  try {
    stan::math::cov_dot_prod(x1, x2, sigma);
  } catch (std::domain_error& e) {
    message = e.what();
  } catch (...) {
    message = "Threw the wrong exection";
  }
  return message;
}


TEST(MathPrimMat, vec_double_cov_dot_prod0) {
  double sigma = 0.5;

  std::vector<double> x(3);
  x[0] = -2;
  x[1] = -1;
  x[2] = -0.5;

  Eigen::MatrixXd cov;
  cov = stan::math::cov_dot_prod(x, sigma);
  EXPECT_NO_THROW(cov = stan::math::cov_dot_prod(x, sigma));
  
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) 
      EXPECT_FLOAT_EQ(sigma + x[i] * x[j], cov(i, j))
        << "index: (" << i << ", " << j << ")";
}
