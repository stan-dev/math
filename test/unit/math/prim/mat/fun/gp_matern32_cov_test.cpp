#include <gtest/gtest.h>
#include <stan/math/prim/mat.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

template <typename T_x, typename T_s, typename T_l>
std::string pull_msg(std::vector<T_x> x, T_s sigma, T_l l) {
  std::string message;
  try {
    stan::math::gp_matern32_cov(x, sigma, l);
  } catch (std::domain_error &e) {
    message = e.what();
  } catch (...) {
    message = "Threw the wrong exception";
  }
  return message;
}

template <typename T_x1, typename T_x2, typename T_s, typename T_l>
std::string pull_msg(std::vector<T_x1> x1, std::vector<T_x2> x2, T_s sigma,
                     T_l l) {
  std::string message;
  try {
    stan::math::gp_matern32_cov(x1, x2, sigma, l);
  } catch (std::domain_error &e) {
    message = e.what();
  } catch (...) {
    message = "Threw the wrong exception";
  }
  return message;
}

TEST(MathPrimMat, vec_double_gp_matern32_cov1) {
  double sigma = 0.2;
  double l = 5.0;

  std::vector<double> x(3);
  x[0] = -2;
  x[1] = -1;
  x[2] = -0.5;

  Eigen::MatrixXd cov;
  cov = stan::math::gp_matern32_cov(x, sigma, l);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      EXPECT_FLOAT_EQ(
          sigma * sigma
              * (1.0
                 + (pow(3.0, 0.5) / l)
                 * stan::math::sqrt(stan::math::squared_distance(x[i], x[j])))
              * std::exp(-1.0 * (pow(3.0, 0.5) / l)
                       * stan::math::sqrt(stan::math::squared_distance(x[i],
                                                                       x[j]))),
          cov(i, j))
          << "index: (" << i << ", " << j << ")";
}

TEST(MathPrimMat, vec_double_ard_gp_matern32_cov1) {
  double sigma = 0.2;
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
  cov = stan::math::gp_matern32_cov(x, sigma, l);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      temp = 0;
      for (int k = 0; k < 5; k++) {
        temp += stan::math::squared_distance(x[i], x[j]) /
          stan::math::square(l[k]);
      }

      EXPECT_FLOAT_EQ(sigma * sigma * (1.0 + pow(3.0, 0.5) * sqrt(temp))
                          * std::exp(-1.0 * pow(3.0, 0.5) * sqrt(temp)),
                      cov(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }
}

TEST(MathPrimMat, vec_eigen_gp_matern32_cov1) {
  double l = 0.2;
  double sigma = 0.3;

  std::vector<Eigen::Matrix<double, 1, -1>> x1(3);
  for (size_t i = 0; i < x1.size(); ++i) {
    x1[i].resize(1, 3);
    x1[i] << 1 * i, 2 * i, 3 * i;
  }

  Eigen::MatrixXd cov;
  cov = stan::math::gp_matern32_cov(x1, sigma, l);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_FLOAT_EQ(
          sigma * sigma
              * (1.0
                 + (pow(3.0, 0.5) / l)
                 * stan::math::sqrt(stan::math::squared_distance(x1[i], x1[j])))
              * std::exp(-1.0 * (pow(3.0, 0.5) / l)
                         * stan::math::sqrt(stan::math::squared_distance(x1[i],
                                                                      x1[j]))),
                         cov(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }
}

TEST(MathPrimMat, vec_eigen_eigen_gp_matern32_cov1) {
  double l = 0.2;
  double sigma = 0.3;

  std::vector<Eigen::Matrix<double, 1, -1>> x1(3);
  for (size_t i = 0; i < x1.size(); ++i) {
    x1[i].resize(1, 3);
    x1[i] << 1 * i, 2 * i, 3 * i;
  }

  std::vector<Eigen::Matrix<double, 1, -1>> x2(3);
  for (size_t i = 0; i < x2.size(); ++i) {
    x2[i].resize(1, 3);
    x2[i] << 2 * i, 3 * i, 4 * i;
  }

  Eigen::MatrixXd cov;
  cov = stan::math::gp_matern32_cov(x1, x2, sigma, l);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_FLOAT_EQ(
          sigma * sigma
              * (1.0
                 + (pow(3.0, 0.5) / l)
                       * stan::math::sqrt(stan::math::squared_distance(x1[i],
                                                                       x2[j])))
          * std::exp(-1.0 * (pow(3.0, 0.5) / l)
                         * stan::math::sqrt(stan::math::squared_distance(x1[i],
                                                                       x2[j]))),
          cov(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }
}

TEST(MathPrimMat, vec_double_double_gp_matern32_cov1) {
  double sigma = 0.2;
  double l = 5.0;

  std::vector<double> x1(3);
  x1[0] = -2;
  x1[1] = -1;
  x1[2] = -0.5;

  std::vector<double> x2(3);
  x2[0] = 3;
  x2[1] = -4;
  x2[2] = -0.7;

  Eigen::MatrixXd cov;
  cov = stan::math::gp_matern32_cov(x1, x2, sigma, l);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      EXPECT_FLOAT_EQ(
          sigma * sigma
              * (1.0
                 + (pow(3.0, 0.5) / l)
                       * stan::math::sqrt(stan::math::squared_distance(x1[i],
                                                                       x2[j])))
              * std::exp(-1.0 * (pow(3.0, 0.5) / l)
                         * stan::math::sqrt(stan::math::squared_distance(x1[i],
                                                                      x2[j]))),
          cov(i, j))
          << "index: (" << i << ", " << j << ")";
}

TEST(MathPrimMat, vec_eigen_ard_gp_matern32_cov1) {
  double sigma = 0.3;
  double temp;

  std::vector<Eigen::Matrix<double, 1, -1>> x1(3);
  for (size_t i = 0; i < x1.size(); ++i) {
    x1[i].resize(1, 3);
    x1[i] << 1 * i, 2 * i, 3 * i;
  }

  std::vector<double> l(5);
  l[0] = 0.1;
  l[1] = 0.2;
  l[2] = 0.3;
  l[3] = 0.4;
  l[4] = 0.5;

  Eigen::MatrixXd cov;
  cov = stan::math::gp_matern32_cov(x1, sigma, l);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      temp = 0;
      for (int k = 0; k < 5; k++) {
        temp += stan::math::squared_distance(x1[i], x1[j]) /
          stan::math::square(l[k]);
      }
      EXPECT_FLOAT_EQ(sigma * sigma * (1.0 + pow(3.0, 0.5) * sqrt(temp))
                          * std::exp(-1.0 * pow(3.0, 0.5) * sqrt(temp)),
                      cov(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }
}

TEST(MathPrimMat, vec_double_double_ard_gp_matern32_cov1) {
  double sigma = 0.2;
  double temp;

  std::vector<double> x1(3);
  x1[0] = -2;
  x1[1] = -1;
  x1[2] = -0.5;

  std::vector<double> x2(3);
  x2[0] = 3;
  x2[1] = -4;
  x2[2] = -0.7;

  std::vector<double> l(4);
  l[0] = 1;
  l[1] = 2;
  l[2] = 3;
  l[3] = 4;

  Eigen::MatrixXd cov;
  cov = stan::math::gp_matern32_cov(x1, x2, sigma, l);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      temp = 0;
      for (int k = 0; k < 4; k++) {
        temp += stan::math::squared_distance(x1[i], x2[j]) /
          stan::math::square(l[k]);
      }
      EXPECT_FLOAT_EQ(sigma * sigma * (1.0 + pow(3.0, 0.5) * pow(temp, 0.5))
                          * std::exp(-1.0 * pow(3.0, 0.5) * pow(temp, 0.5)),
                      cov(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }
}

TEST(MathPrimMat, vec_eigen_eigen_ard_gp_matern32_cov1) {
  double sigma = 0.2;
  double temp;

  std::vector<Eigen::Matrix<double, 1, -1>> x1(3);
  for (size_t i = 0; i < x1.size(); ++i) {
    x1[i].resize(1, 3);
    x1[i] << 1 * i, 2 * i, 3 * i;
  }

  std::vector<Eigen::Matrix<double, 1, -1>> x2(3);
  for (size_t i = 0; i < x2.size(); ++i) {
    x2[i].resize(1, 3);
    x2[i] << 2 * i, 3 * i, 4 * i;
  }

  std::vector<double> l(4);
  l[0] = 1;
  l[1] = 2;
  l[2] = 3;
  l[3] = 4;

  Eigen::MatrixXd cov;
  cov = stan::math::gp_matern32_cov(x1, x2, sigma, l);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      temp = 0;
      for (int k = 0; k < 4; k++) {
        temp += stan::math::squared_distance(x1[i], x2[j]) /
          stan::math::square(l[k]);
      }
      EXPECT_FLOAT_EQ(sigma * sigma * (1.0 + pow(3.0, 0.5) * pow(temp, 0.5))
                          * std::exp(-1.0 * pow(3.0, 0.5) * pow(temp, 0.5)),
                      cov(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }
}

TEST(MathPrimMat, rvec_eigen_gp_matern32_cov1) {
  double l = 0.2;
  double sigma = 0.3;

  std::vector<Eigen::Matrix<double, -1, 1>> x1(3);
  for (size_t i = 0; i < x1.size(); ++i) {
    x1[i].resize(3, 1);
    x1[i] << 1 * i, 2 * i, 3 * i;
  }

  Eigen::MatrixXd cov;
  cov = stan::math::gp_matern32_cov(x1, sigma, l);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_FLOAT_EQ(
          sigma * sigma
              * (1.0
                 + (pow(3.0, 0.5) / l)
                 * stan::math::sqrt(stan::math::squared_distance(x1[i], x1[j])))
              * std::exp(-1.0 * (pow(3.0, 0.5) / l)
                         * stan::math::sqrt(stan::math::squared_distance(x1[i],
                                                                      x1[j]))),
          cov(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }
}

TEST(MathPrimMat, rvec_eigen_ard_gp_matern32_cov1) {
  double sigma = 0.3;
  double temp;

  std::vector<Eigen::Matrix<double, -1, 1>> x1(3);
  for (size_t i = 0; i < x1.size(); ++i) {
    x1[i].resize(3, 1);
    x1[i] << 1 * i, 2 * i, 3 * i;
  }

  std::vector<double> l(4);
  l[0] = 1;
  l[1] = 2;
  l[2] = 3;
  l[3] = 4;

  Eigen::MatrixXd cov;
  cov = stan::math::gp_matern32_cov(x1, sigma, l);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      temp = 0;
      for (int k = 0; k < 4; k++) {
        temp += stan::math::squared_distance(x1[i], x1[j]) /
          stan::math::square(l[k]);
      }

      EXPECT_FLOAT_EQ(sigma * sigma * (1.0 + pow(3.0, 0.5) * sqrt(temp))
                          * std::exp(-1.0 * pow(3.0, 0.5) * sqrt(temp)),
                      cov(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }
}

TEST(MathPrimMat, vec_eigen_rvec_eigen_gp_matern32_cov1) {
  double sigma = 0.3;
  double l = 2.3;

  std::vector<Eigen::Matrix<double, -1, 1>> x1(3);
  for (size_t i = 0; i < x1.size(); ++i) {
    x1[i].resize(3, 1);
    x1[i] << 1 * i, 2 * i, 3 * i;
  }

  std::vector<Eigen::Matrix<double, 1, -1>> x2(3);
  for (size_t i = 0; i < x2.size(); ++i) {
    x2[i].resize(1, 3);
    x2[i] << 2 * i, 3 * i, 4 * i;
  }

  Eigen::MatrixXd cov;
  cov = stan::math::gp_matern32_cov(x1, x2, sigma, l);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_FLOAT_EQ(
          sigma * sigma
              * (1.0
                 + (pow(3.0, 0.5) / l)
                 * stan::math::sqrt(stan::math::squared_distance(x1[i],
                                                                 x2[j])))
              * std::exp(-1.0 * (pow(3.0, 0.5) / l)
                         * stan::math::sqrt(stan::math::squared_distance(x1[i],
                                                                       x2[j]))),
          cov(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }

  cov = stan::math::gp_matern32_cov(x2, x1, sigma, l);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_FLOAT_EQ(
          sigma * sigma
              * (1.0
                 + (pow(3.0, 0.5) / l)
                 * stan::math::sqrt(stan::math::squared_distance(x2[i], x1[j])))
              * std::exp(-1.0 * (pow(3.0, 0.5) / l)
                         * stan::math::sqrt(stan::math::squared_distance(x2[i],
                                                                      x1[j]))),
          cov(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }
}

TEST(MathPrimMat, vec_eigen_rvec_eigen_ard_gp_matern32_cov1) {
  double sigma = 0.3;
  double temp = 0;
  std::vector<double> l(3);
  l[0] = 1;
  l[1] = 2;
  l[2] = 3;

  std::vector<Eigen::Matrix<double, -1, 1>> x1(3);
  for (size_t i = 0; i < x1.size(); ++i) {
    x1[i].resize(3, 1);
    x1[i] << 1 * i, 2 * i, 3 * i;
  }

  std::vector<Eigen::Matrix<double, 1, -1>> x2(3);
  for (size_t i = 0; i < x2.size(); ++i) {
    x2[i].resize(1, 3);
    x2[i] << 2 * i, 3 * i, 4 * i;
  }

  Eigen::MatrixXd cov;
  cov = stan::math::gp_matern32_cov(x1, x2, sigma, l);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      temp = 0;
      for (int k = 0; k < 3; k++) {
        temp += stan::math::squared_distance(x1[i], x2[j]) /
          stan::math::square(l[k]);
      }

      EXPECT_FLOAT_EQ(sigma * sigma * (1.0 + pow(3.0, 0.5) * sqrt(temp))
                          * std::exp(-1.0 * pow(3.0, 0.5) * sqrt(temp)),
                      cov(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }
  cov = stan::math::gp_matern32_cov(x2, x1, sigma, l);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      temp = 0;
      for (int k = 0; k < 3; k++) {
        temp += stan::math::squared_distance(x2[i], x1[j]) /
          stan::math::square(l[k]);
      }

      EXPECT_FLOAT_EQ(sigma * sigma * (1.0 + pow(3.0, 0.5) * sqrt(temp))
                          * std::exp(-1.0 * pow(3.0, 0.5) * sqrt(temp)),
                      cov(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }
}

TEST(MathPrimMat, rvec_eigen_rvec_eigen_gp_matern32_cov1) {
  double sigma = 0.3;
  double l = 2.3;

  std::vector<Eigen::Matrix<double, -1, 1>> x1(3);
  for (size_t i = 0; i < x1.size(); ++i) {
    x1[i].resize(3, 1);
    x1[i] << 1 * i, 2 * i, 3 * i;
  }

  std::vector<Eigen::Matrix<double, -1, 1>> x2(3);
  for (size_t i = 0; i < x2.size(); ++i) {
    x2[i].resize(3, 1);
    x2[i] << 2 * i, 3 * i, 4 * i;
  }

  Eigen::MatrixXd cov;
  cov = stan::math::gp_matern32_cov(x1, x2, sigma, l);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_FLOAT_EQ(
          sigma * sigma
              * (1.0
                 + (pow(3.0, 0.5) / l)
                 * stan::math::sqrt(stan::math::squared_distance(x1[i],
                                                                 x2[j])))
              * std::exp(-1.0 * (pow(3.0, 0.5) / l)
                         * stan::math::sqrt(stan::math::squared_distance(x1[i],
                                                                      x2[j]))),
          cov(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }

  cov = stan::math::gp_matern32_cov(x2, x1, sigma, l);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_FLOAT_EQ(
          sigma * sigma
              * (1.0
                 + (pow(3.0, 0.5) / l)
                 * stan::math::sqrt(stan::math::squared_distance(x2[i],
                                                                 x1[j])))
              * std::exp(-1.0 * (pow(3.0, 0.5) / l)
                         * stan::math::sqrt(stan::math::squared_distance(x2[i],
                                                                      x1[j]))),
          cov(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }
}

TEST(MathPrimMat, rvec_eigen_rvec_eigen_ard_gp_matern32_cov1) {
  double sigma = 0.3;
  double temp;

  std::vector<double> l(3);
  l[0] = 1;
  l[1] = 2;
  l[2] = 3;

  std::vector<Eigen::Matrix<double, -1, 1>> x1(3);
  for (size_t i = 0; i < x1.size(); ++i) {
    x1[i].resize(3, 1);
    x1[i] << 1 * i, 2 * i, 3 * i;
  }

  std::vector<Eigen::Matrix<double, -1, 1>> x2(3);
  for (size_t i = 0; i < x2.size(); ++i) {
    x2[i].resize(3, 1);
    x2[i] << 2 * i, 3 * i, 4 * i;
  }

  Eigen::MatrixXd cov;
  cov = stan::math::gp_matern32_cov(x1, x2, sigma, l);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      temp = 0;
      for (int k = 0; k < 3; k++) {
        temp += stan::math::squared_distance(x1[i], x2[j]) /
          stan::math::square(l[k]);
      }

      EXPECT_FLOAT_EQ(sigma * sigma * (1.0 + pow(3.0, 0.5) * sqrt(temp))
                          * std::exp(-1.0 * pow(3.0, 0.5) * sqrt(temp)),
                      cov(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }
  cov = stan::math::gp_matern32_cov(x2, x1, sigma, l);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      temp = 0;
      for (int k = 0; k < 3; k++) {
        temp += stan::math::squared_distance(x2[i], x1[j]) /
          stan::math::square(l[k]);
      }

      EXPECT_FLOAT_EQ(sigma * sigma * (1.0 + pow(3.0, 0.5) * sqrt(temp))
                          * std::exp(-1.0 * pow(3.0, 0.5) * sqrt(temp)),
                      cov(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }
}

TEST(MathPrimMat, vec_eigen_mixed_gp_matern32_cov2) {
  using stan::math::squared_distance;
  double sigma = 0.2;
  double l = 5;

  std::vector<Eigen::Matrix<double, 1, -1>> x1_rvec(3);
  std::vector<Eigen::Matrix<double, 1, -1>> x2_rvec(4);

  for (size_t i = 0; i < x1_rvec.size(); ++i) {
    x1_rvec[i].resize(1, 3);
    x1_rvec[i] << 1 * i, 2 * i, 3 * i;
  }

  for (size_t i = 0; i < x2_rvec.size(); ++i) {
    x2_rvec[i].resize(1, 3);
    x2_rvec[i] << 2 * i, 3 * i, 4 * i;
  }

  std::vector<Eigen::Matrix<double, -1, 1>> x1_vec(3);
  std::vector<Eigen::Matrix<double, -1, 1>> x2_vec(4);

  for (size_t i = 0; i < x1_vec.size(); ++i) {
    x1_vec[i].resize(3, 1);
    x1_vec[i] << 1 * i, 2 * i, 3 * i;
  }

  for (size_t i = 0; i < x2_vec.size(); ++i) {
    x2_vec[i].resize(3, 1);
    x2_vec[i] << 2 * i, 3 * i, 4 * i;
  }
  Eigen::MatrixXd cov;
  EXPECT_NO_THROW(cov
                  = stan::math::gp_matern32_cov(x1_rvec, x2_vec, sigma, l));
  EXPECT_EQ(3, cov.rows());
  EXPECT_EQ(4, cov.cols());
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_FLOAT_EQ(
          sigma * sigma
              * (1.0
                 + (pow(3.0, 0.5) / l)
                 * stan::math::sqrt(stan::math::squared_distance(x1_rvec[i],
                                                                 x2_vec[j])))
              * std::exp(-1.0 * (pow(3.0, 0.5) / l)
              * stan::math::sqrt(stan::math::squared_distance(x1_rvec[i],
                                                              x2_vec[j]))),
          cov(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }

  Eigen::MatrixXd cov7;
  EXPECT_NO_THROW(cov7
                  = stan::math::gp_matern32_cov(x2_vec, x1_rvec, sigma, l));
  EXPECT_EQ(4, cov7.rows());
  EXPECT_EQ(3, cov7.cols());
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_FLOAT_EQ(
          sigma * sigma
              * (1.0
                 + (pow(3.0, 0.5) / l)
                 * stan::math::sqrt(stan::math::squared_distance(x2_vec[i],
                                                                 x1_rvec[j])))
              * std::exp(-1.0 * (pow(3.0, 0.5) / l)
                      * stan::math::sqrt(stan::math::squared_distance(x2_vec[i],
                                                                  x1_rvec[j]))),
          cov7(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }

  Eigen::MatrixXd cov2;
  EXPECT_NO_THROW(cov2
                  = stan::math::gp_matern32_cov(x1_vec, x2_rvec, sigma, l));
  EXPECT_EQ(3, cov2.rows());
  EXPECT_EQ(4, cov2.cols());
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_FLOAT_EQ(
          sigma * sigma
              * (1.0
                 + (pow(3.0, 0.5) / l)
                     * stan::math::sqrt(stan::math::squared_distance(x1_vec[i],
                                                                  x2_rvec[j])))
              * std::exp(-1.0 * (pow(3.0, 0.5) / l)
                      * stan::math::sqrt(stan::math::squared_distance(x1_vec[i],
                                                                  x2_rvec[j]))),
          cov2(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }

  Eigen::MatrixXd cov8;
  EXPECT_NO_THROW(cov8
                  = stan::math::gp_matern32_cov(x2_rvec, x1_vec, sigma, l));
  EXPECT_EQ(4, cov8.rows());
  EXPECT_EQ(3, cov8.cols());
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_FLOAT_EQ(
          sigma * sigma
              * (1.0
                 + (pow(3.0, 0.5) / l)
                 * stan::math::sqrt(stan::math::squared_distance(x2_rvec[i],
                                                                 x1_vec[j])))
              * std::exp(-1.0 * (pow(3.0, 0.5) / l)
                    * stan::math::sqrt(stan::math::squared_distance(x2_rvec[i],
                                                                   x1_vec[j]))),
          cov8(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }

  Eigen::MatrixXd cov3;
  EXPECT_NO_THROW(cov3
                  = stan::math::gp_matern32_cov(x2_vec, x2_rvec, sigma, l));
  EXPECT_EQ(4, cov3.rows());
  EXPECT_EQ(4, cov3.cols());
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_FLOAT_EQ(
          sigma * sigma
              * (1.0
                 + (pow(3.0, 0.5) / l)
                 * stan::math::sqrt(stan::math::squared_distance(x2_vec[i],
                                                                 x2_rvec[j])))
              * std::exp(-1.0 * (pow(3.0, 0.5) / l)
                      * stan::math::sqrt(stan::math::squared_distance(x2_vec[i],
                                                                  x2_rvec[j]))),
          cov3(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }

  Eigen::MatrixXd cov4;
  EXPECT_NO_THROW(cov4
                  = stan::math::gp_matern32_cov(x2_rvec, x2_vec, sigma, l));
  EXPECT_EQ(4, cov4.rows());
  EXPECT_EQ(4, cov4.cols());
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_FLOAT_EQ(
          sigma * sigma
              * (1.0
                 + (pow(3.0, 0.5) / l)
                 * stan::math::sqrt(stan::math::squared_distance(x2_rvec[i],
                                                                 x2_vec[j])))
              * std::exp(-1.0 * (pow(3.0, 0.5) / l)
                   * stan::math::sqrt(stan::math::squared_distance(x2_rvec[i],
                                                                   x2_vec[j]))),
          cov4(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }

  Eigen::MatrixXd cov5;
  EXPECT_NO_THROW(cov5
                  = stan::math::gp_matern32_cov(x1_rvec, x1_vec, sigma, l));
  EXPECT_EQ(3, cov5.rows());
  EXPECT_EQ(3, cov5.cols());
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_FLOAT_EQ(
          sigma * sigma
              * (1.0
                 + (pow(3.0, 0.5) / l)
                 * stan::math::sqrt(stan::math::squared_distance(x1_rvec[i],
                                                                 x1_vec[j])))
              * std::exp(-1.0 * (pow(3.0, 0.5) / l)
                  * stan::math::sqrt(stan::math::squared_distance(x1_rvec[i],
                                                                  x1_vec[j]))),
          cov5(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }

  Eigen::MatrixXd cov6;
  EXPECT_NO_THROW(cov6
                  = stan::math::gp_matern32_cov(x1_vec, x1_rvec, sigma, l));
  EXPECT_EQ(3, cov6.rows());
  EXPECT_EQ(3, cov6.cols());
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_FLOAT_EQ(
          sigma * sigma
              * (1.0
                 + (pow(3.0, 0.5) / l)
                 * stan::math::sqrt(stan::math::squared_distance(x1_vec[i],
                                                                 x1_rvec[j])))
              * std::exp(-1.0 * (pow(3.0, 0.5) / l)
                      * stan::math::sqrt(stan::math::squared_distance(x1_vec[i],
                                                                  x1_rvec[j]))),
          cov6(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }
}

TEST(MathPrimMat, vec_eigen_mixed_ard_gp_matern32_cov2) {
  using stan::math::squared_distance;
  double sigma = 0.2;
  double temp;

  std::vector<double> l(3);
  l[0] = 1;
  l[1] = 2;
  l[2] = 3;

  std::vector<Eigen::Matrix<double, 1, -1>> x1_rvec(3);
  std::vector<Eigen::Matrix<double, 1, -1>> x2_rvec(4);

  for (size_t i = 0; i < x1_rvec.size(); ++i) {
    x1_rvec[i].resize(1, 3);
    x1_rvec[i] << 1 * i, 2 * i, 3 * i;
  }

  for (size_t i = 0; i < x2_rvec.size(); ++i) {
    x2_rvec[i].resize(1, 3);
    x2_rvec[i] << 2 * i, 3 * i, 4 * i;
  }

  std::vector<Eigen::Matrix<double, -1, 1>> x1_vec(3);
  std::vector<Eigen::Matrix<double, -1, 1>> x2_vec(4);

  for (size_t i = 0; i < x1_vec.size(); ++i) {
    x1_vec[i].resize(3, 1);
    x1_vec[i] << 1 * i, 2 * i, 3 * i;
  }

  for (size_t i = 0; i < x2_vec.size(); ++i) {
    x2_vec[i].resize(3, 1);
    x2_vec[i] << 2 * i, 3 * i, 4 * i;
  }
  Eigen::MatrixXd cov;
  EXPECT_NO_THROW(cov
                  = stan::math::gp_matern32_cov(x1_rvec, x2_vec, sigma, l));
  EXPECT_EQ(3, cov.rows());
  EXPECT_EQ(4, cov.cols());
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      temp = 0;
      for (int k = 0; k < 3; k++) {
        temp += stan::math::squared_distance(x1_rvec[i], x2_vec[j]) /
          stan::math::square(l[k]);
      }

      EXPECT_FLOAT_EQ(sigma * sigma * (1.0 + pow(3.0, 0.5) * sqrt(temp))
                          * std::exp(-1.0 * pow(3.0, 0.5) * sqrt(temp)),
                      cov(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }

  Eigen::MatrixXd cov7;
  EXPECT_NO_THROW(cov7
                  = stan::math::gp_matern32_cov(x2_vec, x1_rvec, sigma, l));
  EXPECT_EQ(4, cov7.rows());
  EXPECT_EQ(3, cov7.cols());
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      temp = 0;
      for (int k = 0; k < 3; k++) {
        temp += stan::math::squared_distance(x2_vec[i], x1_rvec[j]) /
          stan::math::square(l[k]);
      }

      EXPECT_FLOAT_EQ(sigma * sigma * (1.0 + pow(3.0, 0.5) * sqrt(temp))
                          * std::exp(-1.0 * pow(3.0, 0.5) * sqrt(temp)),
                      cov7(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }

  Eigen::MatrixXd cov2;
  EXPECT_NO_THROW(cov2
                  = stan::math::gp_matern32_cov(x1_vec, x2_rvec, sigma, l));
  EXPECT_EQ(3, cov2.rows());
  EXPECT_EQ(4, cov2.cols());
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      temp = 0;
      for (int k = 0; k < 3; k++) {
        temp += stan::math::squared_distance(x1_vec[i], x2_rvec[j]) /
          stan::math::square(l[k]);
      }

      EXPECT_FLOAT_EQ(sigma * sigma * (1.0 + pow(3.0, 0.5) * sqrt(temp))
                          * std::exp(-1.0 * pow(3.0, 0.5) * sqrt(temp)),
                      cov2(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }

  Eigen::MatrixXd cov8;
  EXPECT_NO_THROW(cov8
                  = stan::math::gp_matern32_cov(x2_rvec, x1_vec, sigma, l));
  EXPECT_EQ(4, cov8.rows());
  EXPECT_EQ(3, cov8.cols());
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      temp = 0;
      for (int k = 0; k < 3; k++) {
        temp += stan::math::squared_distance(x2_rvec[i], x1_vec[j]) /
          stan::math::square(l[k]);
      }

      EXPECT_FLOAT_EQ(sigma * sigma * (1.0 + pow(3.0, 0.5) * sqrt(temp))
                          * std::exp(-1.0 * pow(3.0, 0.5) * sqrt(temp)),
                      cov8(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }

  Eigen::MatrixXd cov3;
  EXPECT_NO_THROW(cov3
                  = stan::math::gp_matern32_cov(x2_vec, x2_rvec, sigma, l));
  EXPECT_EQ(4, cov3.rows());
  EXPECT_EQ(4, cov3.cols());
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      temp = 0;
      for (int k = 0; k < 3; k++) {
        temp += stan::math::squared_distance(x2_vec[i], x2_rvec[j]) /
          stan::math::square(l[k]);
      }

      EXPECT_FLOAT_EQ(sigma * sigma * (1.0 + pow(3.0, 0.5) * sqrt(temp))
                          * std::exp(-1.0 * pow(3.0, 0.5) * sqrt(temp)),
                      cov3(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }

  Eigen::MatrixXd cov4;
  EXPECT_NO_THROW(cov4
                  = stan::math::gp_matern32_cov(x2_rvec, x2_vec, sigma, l));
  EXPECT_EQ(4, cov4.rows());
  EXPECT_EQ(4, cov4.cols());
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      temp = 0;
      for (int k = 0; k < 3; k++) {
        temp += stan::math::squared_distance(x2_rvec[i], x2_vec[j]) /
          stan::math::square(l[k]);
      }

      EXPECT_FLOAT_EQ(sigma * sigma * (1.0 + pow(3.0, 0.5) * sqrt(temp))
                          * std::exp(-1.0 * pow(3.0, 0.5) * sqrt(temp)),
                      cov4(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }

  Eigen::MatrixXd cov5;
  EXPECT_NO_THROW(cov5
                  = stan::math::gp_matern32_cov(x1_rvec, x1_vec, sigma, l));
  EXPECT_EQ(3, cov5.rows());
  EXPECT_EQ(3, cov5.cols());
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      temp = 0;
      for (int k = 0; k < 3; k++) {
        temp += stan::math::squared_distance(x1_rvec[i], x1_vec[j]) /
          stan::math::square(l[k]);
      }

      EXPECT_FLOAT_EQ(sigma * sigma * (1.0 + pow(3.0, 0.5) * sqrt(temp))
                          * std::exp(-1.0 * pow(3.0, 0.5) * sqrt(temp)),
                      cov5(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }

  Eigen::MatrixXd cov6;
  EXPECT_NO_THROW(cov6
                  = stan::math::gp_matern32_cov(x1_vec, x1_rvec, sigma, l));
  EXPECT_EQ(3, cov6.rows());
  EXPECT_EQ(3, cov6.cols());
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      temp = 0;
      for (int k = 0; k < 3; k++) {
        temp += stan::math::squared_distance(x1_vec[i], x1_rvec[j]) /
          stan::math::square(l[k]);
      }

      EXPECT_FLOAT_EQ(sigma * sigma * (1.0 + pow(3.0, 0.5) * sqrt(temp))
                          * std::exp(-1.0 * pow(3.0, 0.5) * sqrt(temp)),
                      cov6(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }
}

TEST(MathPrimMat, domain_err_training_sig_l_gamma) {
  double sigma = .2;
  double l = 7;
  double sigma_bad = -1.0;
  double l_bad = -1.0;

  std::vector<double> l_vec(2);
  l_vec[0] = 1;
  l_vec[1] = 2;

  std::vector<double> l_vec_bad(2);
  l_vec_bad[0] = 1;
  l_vec_bad[1] = -1;

  std::vector<double> x(3);
  x[0] = -2;
  x[1] = -1;
  x[2] = -0.5;

  std::vector<Eigen::Matrix<double, -1, 1>> x_2(3);
  for (size_t i = 0; i < x_2.size(); ++i) {
    x_2[i].resize(3, 1);
    x_2[i] << 1, 2, 3;
  }

  std::vector<Eigen::Matrix<double, 1, -1>> x_3(3);
  for (size_t i = 0; i < x_3.size(); ++i) {
    x_3[i].resize(1, 3);
    x_3[i] << 1, 2, 3;
  }

  std::string msg1, msg2, msg3, msg4;
  msg1 = pull_msg(x, sigma, l_bad);
  msg2 = pull_msg(x, sigma_bad, l_bad);
  msg3 = pull_msg(x, sigma_bad, l_bad);
  EXPECT_TRUE(std::string::npos != msg1.find(" length scale")) << msg1;
  EXPECT_TRUE(std::string::npos != msg2.find(" marginal variance")) << msg2;
  EXPECT_TRUE(std::string::npos != msg3.find(" marginal variance")) << msg3;

  msg1 = pull_msg(x, sigma, l_bad);
  msg2 = pull_msg(x, sigma_bad, l_vec);
  msg3 = pull_msg(x, sigma_bad, l_bad);
  EXPECT_TRUE(std::string::npos != msg1.find(" length scale")) << msg1;
  EXPECT_TRUE(std::string::npos != msg2.find(" marginal variance")) << msg2;
  EXPECT_TRUE(std::string::npos != msg3.find(" marginal variance")) << msg3;

  msg1 = pull_msg(x, sigma, l_vec_bad);
  msg2 = pull_msg(x, sigma_bad, l_vec_bad);
  msg3 = pull_msg(x, sigma_bad, l_vec);
  EXPECT_TRUE(std::string::npos != msg1.find(" length scale")) << msg1;
  EXPECT_TRUE(std::string::npos != msg2.find(" marginal variance")) << msg2;
  EXPECT_TRUE(std::string::npos != msg3.find(" marginal variance")) << msg3;

  msg1 = pull_msg(x_2, sigma, l_bad);
  msg2 = pull_msg(x_2, sigma_bad, l);
  msg3 = pull_msg(x_2, sigma_bad, l_bad);
  EXPECT_TRUE(std::string::npos != msg1.find(" length scale")) << msg1;
  EXPECT_TRUE(std::string::npos != msg2.find(" marginal variance")) << msg2;
  EXPECT_TRUE(std::string::npos != msg3.find(" marginal variance")) << msg3;

  msg1 = pull_msg(x_2, sigma, l_bad);
  msg2 = pull_msg(x_2, sigma_bad, l_vec);
  msg3 = pull_msg(x_2, sigma_bad, l_bad);
  EXPECT_TRUE(std::string::npos != msg1.find(" length scale")) << msg1;
  EXPECT_TRUE(std::string::npos != msg2.find(" marginal variance")) << msg2;
  EXPECT_TRUE(std::string::npos != msg3.find(" marginal variance")) << msg3;

  msg1 = pull_msg(x_2, sigma, l_vec_bad);
  msg2 = pull_msg(x_2, sigma_bad, l_vec_bad);
  msg3 = pull_msg(x_2, sigma_bad, l_vec);
  EXPECT_TRUE(std::string::npos != msg1.find(" length scale")) << msg1;
  EXPECT_TRUE(std::string::npos != msg2.find(" marginal variance")) << msg2;
  EXPECT_TRUE(std::string::npos != msg3.find(" marginal variance")) << msg3;

  msg1 = pull_msg(x_2, x_3, sigma, l_bad);
  msg2 = pull_msg(x_2, x_3, sigma_bad, l);
  msg3 = pull_msg(x_2, x_3, sigma_bad, l_bad);
  EXPECT_TRUE(std::string::npos != msg1.find(" length scale")) << msg1;
  EXPECT_TRUE(std::string::npos != msg2.find(" marginal variance")) << msg2;
  EXPECT_TRUE(std::string::npos != msg3.find(" marginal variance")) << msg3;

  msg1 = pull_msg(x_2, x_3, sigma, l_bad);
  msg2 = pull_msg(x_2, x_3, sigma_bad, l_vec);
  msg3 = pull_msg(x_2, x_3, sigma_bad, l_bad);
  EXPECT_TRUE(std::string::npos != msg1.find(" length scale")) << msg1;
  EXPECT_TRUE(std::string::npos != msg2.find(" marginal variance")) << msg2;
  EXPECT_TRUE(std::string::npos != msg3.find(" marginal variance")) << msg3;

  msg1 = pull_msg(x_2, x_3, sigma, l_vec_bad);
  msg2 = pull_msg(x_2, x_3, sigma_bad, l_vec_bad);
  msg3 = pull_msg(x_2, x_3, sigma_bad, l_vec);
  EXPECT_TRUE(std::string::npos != msg1.find(" length scale")) << msg1;
  EXPECT_TRUE(std::string::npos != msg2.find(" marginal variance")) << msg2;
  EXPECT_TRUE(std::string::npos != msg3.find(" marginal variance")) << msg3;

  msg1 = pull_msg(x_3, x_3, sigma, l_vec_bad);
  msg2 = pull_msg(x_3, x_3, sigma_bad, l_vec_bad);
  msg3 = pull_msg(x_3, x_3, sigma_bad, l_vec);
  EXPECT_TRUE(std::string::npos != msg1.find(" length scale")) << msg1;
  EXPECT_TRUE(std::string::npos != msg2.find(" marginal variance")) << msg2;
  EXPECT_TRUE(std::string::npos != msg3.find(" marginal variance")) << msg3;

  msg1 = pull_msg(x_2, x_2, sigma, l_vec_bad);
  msg2 = pull_msg(x_2, x_2, sigma_bad, l_vec_bad);
  msg3 = pull_msg(x_2, x_2, sigma_bad, l_vec);
  EXPECT_TRUE(std::string::npos != msg1.find(" length scale")) << msg1;
  EXPECT_TRUE(std::string::npos != msg2.find(" marginal variance")) << msg2;
  EXPECT_TRUE(std::string::npos != msg3.find(" marginal variance")) << msg3;
}

TEST(MathPrimMat, nan_error_training_sig_l_gamma) {
  double sigma = 0.2;
  double l = 5;

  std::vector<double> x(3);
  x[0] = -2;
  x[1] = -1;
  x[2] = -0.5;

  std::vector<Eigen::Matrix<double, -1, 1>> x_2(3);
  for (size_t i = 0; i < x_2.size(); ++i) {
    x_2[i].resize(3, 1);
    x_2[i] << 1, 2, 3;
  }

  std::vector<Eigen::Matrix<double, 1, -1>> x_3(3);
  for (size_t i = 0; i < x_3.size(); ++i) {
    x_3[i].resize(1, 3);
    x_3[i] << 1, 2, 3;
  }

  std::vector<double> l_vec(2);
  l_vec[0] = 1;
  l_vec[1] = 2;

  std::vector<double> l_vec_bad(2);
  l_vec_bad[0] = 1;
  l_vec_bad[1] = -1;

  std::vector<double> x_bad(x);
  x_bad[1] = std::numeric_limits<double>::quiet_NaN();

  std::vector<Eigen::Matrix<double, -1, 1>> x_bad_2(x_2);
  x_bad_2[1](1) = std::numeric_limits<double>::quiet_NaN();

  std::vector<Eigen::Matrix<double, 1, -1>> x_bad_3(x_3);
  x_bad_3[1](1) = std::numeric_limits<double>::quiet_NaN();

  double sigma_bad = std::numeric_limits<double>::quiet_NaN();
  double l_bad = std::numeric_limits<double>::quiet_NaN();

  std::string msg1, msg2, msg3, msg4;
  msg1 = pull_msg(x, sigma, l_bad);
  msg2 = pull_msg(x, sigma_bad, l);
  msg3 = pull_msg(x, sigma_bad, l_bad);
  EXPECT_TRUE(std::string::npos != msg1.find(" length scale")) << msg1;
  EXPECT_TRUE(std::string::npos != msg2.find(" marginal variance")) << msg2;
  EXPECT_TRUE(std::string::npos != msg3.find(" marginal variance")) << msg3;

  EXPECT_THROW(stan::math::gp_matern32_cov(x, sigma, l_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_matern32_cov(x, sigma_bad, l),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_matern32_cov(x, sigma_bad, l_bad),
               std::domain_error);

  EXPECT_THROW(stan::math::gp_matern32_cov(x_bad, l, sigma),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_matern32_cov(x_bad, sigma_bad, l),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_matern32_cov(x_bad, sigma_bad, l_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_matern32_cov(x_bad, sigma, l_bad),
               std::domain_error);

  EXPECT_THROW(stan::math::gp_matern32_cov(x_2, sigma, l_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_matern32_cov(x_2, sigma_bad, l),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_matern32_cov(x_2, sigma_bad, l_bad),
               std::domain_error);

  EXPECT_THROW(stan::math::gp_matern32_cov(x_bad_2, sigma, l_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_matern32_cov(x_bad_2, sigma_bad, l),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_matern32_cov(x_bad_2, sigma, l),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_matern32_cov(x_bad_2, sigma_bad, l_bad),
               std::domain_error);

  EXPECT_THROW(stan::math::gp_matern32_cov(x_3, sigma, l_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_matern32_cov(x_3, sigma_bad, l),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_matern32_cov(x_3, sigma_bad, l_bad),
               std::domain_error);

  EXPECT_THROW(stan::math::gp_matern32_cov(x_bad_3, sigma, l_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_matern32_cov(x_bad_3, sigma_bad, l),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_matern32_cov(x_bad_3, sigma_bad, l_bad),
               std::domain_error);

  EXPECT_TRUE(std::string::npos != msg1.find(" length scale")) << msg1;
  EXPECT_TRUE(std::string::npos != msg2.find(" marginal variance")) << msg2;
  EXPECT_TRUE(std::string::npos != msg3.find(" marginal variance")) << msg3;

  EXPECT_THROW(stan::math::gp_matern32_cov(x, sigma, l_vec_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_matern32_cov(x, sigma_bad, l_vec),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_matern32_cov(x, sigma_bad, l_vec_bad),
               std::domain_error);

  EXPECT_THROW(stan::math::gp_matern32_cov(x_bad, sigma, l_vec),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_matern32_cov(x_bad, sigma_bad, l_vec),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_matern32_cov(x_bad, sigma_bad, l_vec_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_matern32_cov(x_bad, sigma, l_vec_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_matern32_cov(x_bad, sigma_bad, l_vec_bad),
               std::domain_error);

  EXPECT_THROW(stan::math::gp_matern32_cov(x_2, sigma, l_vec_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_matern32_cov(x_2, sigma_bad, l_vec),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_matern32_cov(x_2, sigma_bad, l_vec_bad),
               std::domain_error);

  EXPECT_THROW(stan::math::gp_matern32_cov(x_bad_2, sigma, l_vec_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_matern32_cov(x_bad_2, sigma_bad, l_vec),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_matern32_cov(x_bad_2, sigma_bad, l_vec_bad),
               std::domain_error);

  EXPECT_THROW(stan::math::gp_matern32_cov(x_3, sigma, l_vec_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_matern32_cov(x_3, sigma_bad, l_vec),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_matern32_cov(x_3, sigma_bad, l_vec_bad),
               std::domain_error);

  EXPECT_THROW(stan::math::gp_matern32_cov(x_bad_3, sigma, l_vec_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_matern32_cov(x_bad_3, sigma_bad, l_vec),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_matern32_cov(x_bad_3, sigma_bad, l_vec_bad),
               std::domain_error);
}

TEST(MathPrimMat, dim_mismatch_vec_eigen_rvec_gp_matern32_cov2) {
  double sigma = 0.2;
  double l = 5;

  std::vector<Eigen::Matrix<double, 1, -1>> x_vec_1(3);
  for (size_t i = 0; i < x_vec_1.size(); ++i) {
    x_vec_1[i].resize(1, 3);
    x_vec_1[i] << 1, 2, 3;
  }

  std::vector<Eigen::Matrix<double, 1, -1>> x_vec_2(4);
  for (size_t i = 0; i < x_vec_2.size(); ++i) {
    x_vec_2[i].resize(1, 4);
    x_vec_2[i] << 4, 1, 3, 1;
  }
  EXPECT_THROW(stan::math::gp_matern32_cov(x_vec_1, x_vec_2, sigma, l),
               std::invalid_argument);
}

TEST(MathPrimMat, dim_mismatch_vec_eigen_mixed_gp_matern32_cov) {
  double sigma = 0.2;
  double l = 5;

  std::vector<Eigen::Matrix<double, 1, -1>> x_rvec_1(3);
  for (size_t i = 0; i < x_rvec_1.size(); ++i) {
    x_rvec_1[i].resize(1, 3);
    x_rvec_1[i] << 1, 2, 3;
  }

  std::vector<Eigen::Matrix<double, -1, 1>> x_vec_1(4);
  for (size_t i = 0; i < x_vec_1.size(); ++i) {
    x_vec_1[i].resize(4, 1);
    x_vec_1[i] << 4, 1, 3, 1;
  }

  std::vector<Eigen::Matrix<double, -1, 1>> x_vec_2(3);
  for (size_t i = 0; i < x_vec_2.size(); ++i) {
    x_vec_2[i].resize(3, 1);
    x_vec_2[i] << 1, 2, 3;
  }

  std::vector<Eigen::Matrix<double, 1, -1>> x_rvec_2(4);
  for (size_t i = 0; i < x_rvec_2.size(); ++i) {
    x_rvec_2[i].resize(1, 4);
    x_rvec_2[i] << 1, 2, 3, 4;
  }
  EXPECT_THROW(stan::math::gp_matern32_cov(x_rvec_1, x_vec_1, sigma, l),
               std::invalid_argument);
  EXPECT_THROW(stan::math::gp_matern32_cov(x_vec_1, x_rvec_1, sigma, l),
               std::invalid_argument);
  EXPECT_THROW(stan::math::gp_matern32_cov(x_rvec_2, x_vec_2, sigma, l),
               std::invalid_argument);
  EXPECT_THROW(stan::math::gp_matern32_cov(x_vec_2, x_rvec_2, sigma, l),
               std::invalid_argument);
}
