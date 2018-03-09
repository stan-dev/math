#include <gtest/gtest.h>
#include <limits>
#include <stan/math/prim/mat.hpp>
#include <string>
#include <vector>

template <typename T_x, typename T_sigma>
std::string pull_msg(std::vector<T_x> x, T_sigma sigma) {
  std::string message;
  try {
    stan::math::gp_dot_prod_cov(x, sigma);
  } catch (std::domain_error &e) {
    message = e.what();
  } catch (...) {
    message = "Threw the wrong exection";
  }
  return message;
}

template <typename T_x1, typename T_x2, typename T_sigma>
std::string pull_msg(std::vector<T_x1> x1, std::vector<T_x2> x2,
                     T_sigma sigma) {
  std::string message;
  try {
    stan::math::gp_dot_prod_cov(x1, x2, sigma);
  } catch (std::domain_error &e) {
    message = e.what();
  } catch (...) {
    message = "Threw the wrong exection";
  }
  return message;
}

TEST(MathPrimMat, dummy_function) {
  // just testing templates here

  std::vector<double> x_0(2);
  x_0[0] = 1;
  x_0[1] = 2;

  std::vector<Eigen::Matrix<double, -1, 1>> x_1(2);
  for (size_t i = 0; i < x_1.size(); ++i) {
    x_1[i].resize(3, 1);
    x_1[i] << 2 * i, 3 * i, 4 * i;
  }

  // stan::math::dummy_fun(x_0);
  // stan::math::dummy_fun(x_1);
}

TEST(MathPrimMat, vec_double_gp_dot_prod_cov0) {
  double sigma = 0.5;
  double sigma_sq = pow(sigma, 2);

  std::vector<double> x(3);
  x[0] = -2;
  x[1] = -1;
  x[2] = -0.5;

  Eigen::MatrixXd cov;
  cov = stan::math::gp_dot_prod_cov(x, sigma);
  EXPECT_NO_THROW(cov = stan::math::gp_dot_prod_cov(x, sigma));

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      EXPECT_FLOAT_EQ(sigma_sq + x[i] * x[j], cov(i, j))
          << "index: (" << i << ", " << j << ")";
}

TEST(MathPrimMat, vec_x_gp_dot_prod_cov0) {
  double sigma = 0.5;

  std::vector<Eigen::Matrix<double, -1, 1>> x(3);
  for (size_t i = 0; i < x.size(); ++i) {
    x[i].resize(3, 1);
    x[i] << 2 * i, 3 * i, 4 * i;
  }

  Eigen::MatrixXd cov;
  EXPECT_NO_THROW(cov = stan::math::gp_dot_prod_cov(x, sigma));
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_FLOAT_EQ(sigma * sigma + stan::math::dot_product(x[i], x[j]),
                      cov(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }
}

/////////////////////////////////

TEST(MathPrimMat, vec_NaN_x_gp_dot_prod_cov_cov0) {
  double sigma = 1;
  std::vector<double> x(2);
  x[0] = 1;
  x[1] = std::numeric_limits<double>::quiet_NaN();

  EXPECT_THROW(stan::math::gp_dot_prod_cov(x, sigma), std::domain_error);
  std::string msg;
  msg = pull_msg(x, sigma);
  EXPECT_TRUE(std::string::npos != msg.find(" x")) << msg;
}

TEST(MathPrimMat, vec_NaN_sigma_gp_dot_prod_cov_cov0) {
  double sigma = std::numeric_limits<double>::quiet_NaN();
  std::vector<double> x(2);
  x[0] = 1;
  x[1] = 2;

  EXPECT_THROW(stan::math::gp_dot_prod_cov(x, sigma), std::domain_error);
  std::string msg;
  msg = pull_msg(x, sigma);
  EXPECT_TRUE(std::string::npos != msg.find(" sigma")) << msg;
}

TEST(MathPrimMat, vec_NaN_x1_gp_dot_prod_cov_cov0) {
  double sigma = 2;

  std::vector<double> x1(2);
  x1[0] = std::numeric_limits<double>::quiet_NaN();
  x1[1] = 2;

  std::vector<double> x2(2);
  x2[0] = 1;
  x2[1] = 2;

  EXPECT_THROW(stan::math::gp_dot_prod_cov(x1, x2, sigma), std::domain_error);
  std::string msg;
  msg = pull_msg(x1, x2, sigma);
  EXPECT_TRUE(std::string::npos != msg.find(" x1")) << msg;
}

TEST(MathPrimMat, vec_NaN_x2_gp_dot_prod_cov_cov0) {
  double sigma = 2;

  std::vector<double> x1(2);
  x1[0] = 1;
  x1[1] = 2;

  std::vector<double> x2(2);
  x2[0] = 1;
  x2[1] = std::numeric_limits<double>::quiet_NaN();

  EXPECT_THROW(stan::math::gp_dot_prod_cov(x1, x2, sigma), std::domain_error);
  std::string msg;
  msg = pull_msg(x1, x2, sigma);
  EXPECT_TRUE(std::string::npos != msg.find(" x2")) << msg;
}

TEST(MathPrimMat, vec_NaN_x1_x2_sigma_gp_dot_prod_cov_cov0) {
  double sigma = std::numeric_limits<double>::quiet_NaN();

  std::vector<double> x1(2);
  x1[0] = 1;
  x1[1] = 2;

  std::vector<double> x2(2);
  x2[0] = 1;
  x2[1] = 3;

  EXPECT_THROW(stan::math::gp_dot_prod_cov(x1, x2, sigma), std::domain_error);
  std::string msg;
  msg = pull_msg(x1, x2, sigma);
  EXPECT_TRUE(std::string::npos != msg.find(" sigma")) << msg;
}

TEST(MathPrimMat, vec_inf_x_gp_dot_prod_cov_cov0) {
  double sigma = 1;
  std::vector<double> x(2);
  x[0] = 1;
  x[1] = std::numeric_limits<double>::infinity();

  EXPECT_THROW(stan::math::gp_dot_prod_cov(x, sigma), std::domain_error);
  std::string msg;
  msg = pull_msg(x, sigma);
  EXPECT_TRUE(std::string::npos != msg.find(" x")) << msg;
}

TEST(MathPrimMat, vec_inf_sigma_gp_dot_prod_cov_cov0) {
  double sigma = std::numeric_limits<double>::infinity();
  std::vector<double> x(2);
  x[0] = 1;
  x[1] = 2;

  EXPECT_THROW(stan::math::gp_dot_prod_cov(x, sigma), std::domain_error);
  std::string msg;
  msg = pull_msg(x, sigma);
  EXPECT_TRUE(std::string::npos != msg.find(" sigma")) << msg;
}

TEST(MathPrimMat, vec_inf_x1_gp_dot_prod_cov_cov0) {
  double sigma = 2;

  std::vector<double> x1(2);
  x1[0] = std::numeric_limits<double>::infinity();
  x1[1] = 2;

  std::vector<double> x2(2);
  x2[0] = 1;
  x2[1] = 2;

  EXPECT_THROW(stan::math::gp_dot_prod_cov(x1, x2, sigma), std::domain_error);
  std::string msg;
  msg = pull_msg(x1, x2, sigma);
  EXPECT_TRUE(std::string::npos != msg.find(" x1")) << msg;
}

TEST(MathPrimMat, vec_inf_x2_gp_dot_prod_cov_cov0) {
  double sigma = 2;

  std::vector<double> x1(2);
  x1[0] = 1;
  x1[1] = 2;

  std::vector<double> x2(2);
  x2[0] = 1;
  x2[1] = std::numeric_limits<double>::infinity();

  EXPECT_THROW(stan::math::gp_dot_prod_cov(x1, x2, sigma), std::domain_error);
  std::string msg;
  msg = pull_msg(x1, x2, sigma);
  EXPECT_TRUE(std::string::npos != msg.find(" x2")) << msg;
}

TEST(MathPrimMat, vec_inf_x1_x2_sigma_gp_dot_prod_cov_cov0) {
  double sigma = std::numeric_limits<double>::infinity();

  std::vector<double> x1(2);
  x1[0] = 1;
  x1[1] = 2;

  std::vector<double> x2(2);
  x2[0] = 1;
  x2[1] = 3;

  EXPECT_THROW(stan::math::gp_dot_prod_cov(x1, x2, sigma), std::domain_error);

  std::string msg;
  msg = pull_msg(x1, x2, sigma);
  EXPECT_TRUE(std::string::npos != msg.find(" sigma")) << msg;
}

TEST(MathPrimMat, rvec_x_gp_dot_prod_cov0) {
  double sigma = 0.5;

  std::vector<Eigen::Matrix<double, 1, -1>> x(3);
  for (size_t i = 0; i < x.size(); ++i) {
    x[i].resize(1, 3);
    x[i] << 1 * i, 2 * i, 3 * i;
  }

  Eigen::MatrixXd cov;
  EXPECT_NO_THROW(cov = stan::math::gp_dot_prod_cov(x, sigma));
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_FLOAT_EQ(sigma * sigma + stan::math::dot_product(x[i], x[j]),
                      cov(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }
}

TEST(MathPrimMat, vec_vec_x1_x2_gp_dot_prod_cov0) {
  double sigma = 0.5;

  std::vector<Eigen::Matrix<double, -1, 1>> x1(3);
  for (size_t i = 0; i < x1.size(); ++i) {
    x1[i].resize(3, 1);
    x1[i] << 2 * i, 3 * i, 4 * i;
  }

  std::vector<Eigen::Matrix<double, -1, 1>> x2(4);
  for (size_t i = 0; i < x2.size(); ++i) {
    x2[i].resize(3, 1);
    x2[i] << 1 * i, 2 * i, 3 * i;
  }

  Eigen::MatrixXd cov;
  EXPECT_NO_THROW(cov = stan::math::gp_dot_prod_cov(x1, x2, sigma));
  EXPECT_EQ(3, cov.rows());
  EXPECT_EQ(4, cov.cols());
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 4; ++j) {
      EXPECT_FLOAT_EQ(sigma * sigma + stan::math::dot_product(x1[i], x2[j]),
                      cov(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }
}

TEST(MathPrimMat, rvec_vec_x1_x2_gp_dot_prod_cov0) {
  double sigma = 0.5;

  std::vector<Eigen::Matrix<double, 1, -1>> x1(3);
  std::vector<Eigen::Matrix<double, 1, -1>> x2(4);

  for (size_t i = 0; i < x1.size(); ++i) {
    x1[i].resize(1, 3);
    x1[i] << 1 * i, 2 * i, 3 * i;
  }

  for (size_t i = 0; i < x2.size(); ++i) {
    x2[i].resize(1, 3);
    x2[i] << 2 * i, 3 * i, 4 * i;
  }

  Eigen::MatrixXd cov;
  EXPECT_NO_THROW(cov = stan::math::gp_dot_prod_cov(x1, x2, sigma));
  EXPECT_EQ(3, cov.rows());
  EXPECT_EQ(4, cov.cols());
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 4; ++j) {
      EXPECT_FLOAT_EQ(sigma * sigma + stan::math::dot_product(x1[i], x2[j]),
                      cov(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }

  Eigen::MatrixXd cov2;
  EXPECT_NO_THROW(cov2 = stan::math::gp_dot_prod_cov(x1, x2, sigma));
  EXPECT_EQ(3, cov2.rows());
  EXPECT_EQ(4, cov2.cols());
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 4; ++j) {
      EXPECT_FLOAT_EQ(sigma * sigma + stan::math::dot_product(x1[i], x2[j]),
                      cov2(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }
}

TEST(MathPrimMat, vec_rvec_x1_x2_gp_dot_prod_cov0) {
  double sigma = 0.7;

  std::vector<Eigen::Matrix<double, 1, -1>> x1_rvec(3);

  for (size_t i = 0; i < x1_rvec.size(); ++i) {
    x1_rvec[i].resize(1, 3);
    x1_rvec[i] << 1 * i, 2 * i, 3 * i;
  }

  std::vector<Eigen::Matrix<double, -1, 1>> x1_vec(3);

  for (size_t i = 0; i < x1_vec.size(); ++i) {
    x1_vec[i].resize(3, 1);
    x1_vec[i] << 1 * i, 2 * i, 3 * i;
  }

  Eigen::MatrixXd cov;
  EXPECT_NO_THROW(cov = stan::math::gp_dot_prod_cov(x1_vec, x1_rvec, sigma));
  EXPECT_EQ(3, cov.rows());
  EXPECT_EQ(3, cov.cols());
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      EXPECT_FLOAT_EQ(sigma * sigma +
                          stan::math::dot_product(x1_vec[i], x1_rvec[j]),
                      cov(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }
}

TEST(MathPrimMat, rvec_rvec_x1_x2_gp_dot_prod_cov0) {
  double sigma = 0.8;

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

  Eigen::MatrixXd cov;
  EXPECT_NO_THROW(cov = stan::math::gp_dot_prod_cov(x1_rvec, x2_rvec, sigma));
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      EXPECT_FLOAT_EQ(sigma * sigma +
                          stan::math::dot_product(x1_rvec[i], x2_rvec[j]),
                      cov(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }
}
