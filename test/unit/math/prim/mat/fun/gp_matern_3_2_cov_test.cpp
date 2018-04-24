#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <string>
#include <vector>

// temp:
#include <cmath>

template <typename T_x, typename T_l, typename T_s, typename T_g>
std::string pull_msg(std::vector<T_x> x, T_l l, T_s sigma, T_g gamma) {
  std::string message;
  try {
    stan::math::gp_matern_3_2_cov(x, l, sigma, gamma);
  } catch (std::domain_error& e) {
    message = e.what();
  } catch (...) {
    message = "Threw the wrong exception";
  }
  return message;
}

TEST(MathPrimMat, vec_double_gp_matern_3_2_cov1){
  double sigma = 0.2;
  double gamma = 1.2;
  double l = 5.0;

  std::vector<double> x(3);
  x[0] = -2;
  x[1] = -1;
  x[2] = -0.5;

  Eigen::MatrixXd cov;
  cov = stan::math::gp_matern_3_2_cov(x, l, sigma, gamma);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      EXPECT_FLOAT_EQ(sigma * sigma * (1.0 + (pow(3.0, 0.5) / l) * gamma *
                                       abs(x[i] - x[j])) *
                      std::exp(-1.0 * gamma * (pow(3.0, 0.5) / l) *
                               abs(x[i] - x[j])),
                      cov(i, j))
        << "index: (" << i << ", " << j << ")";
}

TEST(MathPrimMat, vec_double_ard_gp_matern_3_2_cov1) {
  double sigma = 0.2;
  double gamma = 1.2;
  double temp; 
  
  std::vector<double> x(3);
  x[0] = -2;
  x[1] = -1;
  x[2] = -0.5;

  std::vector<double> l(5);
  l[0] = 0.1;
  l[1] = 0.2;
  l[2] = 0.3;
  l[3] = 0.4;
  l[4] = 0.5;
  
  Eigen::MatrixXd cov;
  cov = stan::math::gp_matern_3_2_cov(x, l, sigma, gamma);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      temp = 0;
      for (int k = 0; k < 5; k++) {
        temp += (1.0 + gamma * pow(3.0, 0.5) * abs(x[i] - x[j]) / l[k]) *
          exp(-1.0 * gamma * pow(3.0, 0.5) * abs(x[i] - x[j]) / l[k]);
      }
      EXPECT_FLOAT_EQ(sigma * sigma * temp,
                      cov(i, j))
        << "index: (" << i << ", " << j << ")";
    }
  }
}

TEST(MathPrimMat, vec_eigen_gp_matern_3_2_cov1) {
  double l = 0.2;
  double sigma = 0.3;
  double gamma = 0.5;
  double distance = 0;

  std::vector<Eigen::Matrix<double, 1, -1> > x1(3);
  std::vector<Eigen::Matrix<double, 1, -1> > x2(4);
  for (size_t i = 0; i < x1.size(); ++i) {
    x1[i].resize(1, 3);
    x1[i] << 1 * i, 2 * i, 3 * i;
  }

  for (size_t i = 0; i < x2.size(); ++i) {
    x2[i].resize(1, 3);
    x2[i] << 2 * i, 3 * i, 4 * i;
  }
  
  Eigen::MatrixXd cov;
  cov = stan::math::gp_matern_3_2_cov(x1, x2, l, sigma, gamma);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      distance = 0;
      for (int k = 0; k < 3; k++)
        distance += abs(x1[i][k] - x2[j][k]);
      EXPECT_FLOAT_EQ(sigma * sigma * (1.0 + (pow(3.0, 0.5) / l) * gamma *
                                       distance) *
                      std::exp(-1.0 * gamma * (pow(3.0, 0.5) / l) *
                               distance),
                      cov(i, j))
        << "index: (" << i << ", " << j << ")";
    }
  }
}
