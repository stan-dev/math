#include <gtest/gtest.h>
#include <stan/math/prim/mat.hpp>
#include <limits>
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

TEST(MathPrimMat, vec_double_gp_dot_prod_cov0) {
  double sigma = 0.5;

  std::vector<double> x(3);
  x[0] = -2;
  x[1] = -1;
  x[2] = -0.5;

  Eigen::MatrixXd ref(3, 3);

  ref << 1.00, 0.500, 0.2500,
    0.50, 0.250, 0.1250,
    0.25, 0.125, 0.0625;

  Eigen::MatrixXd cov;
  EXPECT_NO_THROW(cov = stan::math::gp_dot_prod_cov(x, sigma));

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      EXPECT_FLOAT_EQ(ref(i, j), cov(i, j))
          << "index: (" << i << ", " << j << ")";
}

TEST(MathPrimMat, vec_vec_double_gp_dot_prod_cov0) {
  double sigma = 0.5;

  std::vector<double> x1(3);
  x1[0] = -2;
  x1[1] = -1;
  x1[2] = -0.5;

  std::vector<double> x2(4);
  x2[0] = 1;
  x2[1] = 2;
  x2[2] = 3;
  x2[3] = 4;

  Eigen::MatrixXd ref(3, 4);
  ref << -0.500, -1.00, -1.500, -2.0,
    -0.250, -0.50, -0.750, -1.0,
    -0.125, -0.25, -0.375, -0.5;

  Eigen::MatrixXd cov;
  EXPECT_NO_THROW(cov = stan::math::gp_dot_prod_cov(x1, x2, sigma));

  EXPECT_EQ(ref.rows(), cov.rows());
  EXPECT_EQ(ref.cols(), cov.cols());
  for (int i = 0; i < ref.rows(); i++)
    for (int j = 0; j < ref.cols(); j++)
      EXPECT_FLOAT_EQ(ref(i, j), cov(i, j))
          << "index: (" << i << ", " << j << ")";
}

TEST(MathPrimMat, vec_x_gp_dot_prod_cov0_bad_size_sigma) {
  Eigen::MatrixXd Sigma(2, 2);

  Sigma << 6.8, 11.2, 11.2, 19.8;

  std::vector<Eigen::Matrix<double, -1, 1>> x(3);
  for (size_t i = 0; i < x.size(); ++i) {
    x[i].resize(3, 1);
    x[i] << 2 * (i + 1), 3 * (i + 1), 4 * (i + 1);
  }

  EXPECT_THROW(stan::math::gp_dot_prod_cov(x, Sigma), std::invalid_argument);
}

TEST(MathPrimMat, vec_x_gp_dot_prod_cov0) {
  Eigen::MatrixXd Sigma(2, 2);

  Sigma << 1.5, 0.5,
    0.5, 1.7;

  Eigen::MatrixXd ref(4, 4);

  ref << 27.3, 54.6, 81.9, 109.2,
    54.6, 109.2, 163.8, 218.4,
    81.9, 163.8, 245.7, 327.6,
    109.2, 218.4, 327.6, 436.8;

  std::vector<Eigen::Matrix<double, -1, 1>> x(4);
  for (size_t i = 0; i < x.size(); ++i) {
    x[i].resize(2, 1);
    x[i] << 2 * (i + 1), 3 * (i + 1);
  }

  Eigen::MatrixXd cov;
  EXPECT_NO_THROW(cov = stan::math::gp_dot_prod_cov(x, Sigma));
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_FLOAT_EQ(ref(i, j), cov(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }
}

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

TEST(MathPrimMat, vec_vec_x1_x2_gp_dot_prod_cov0_bad_size_vectors) {
  Eigen::MatrixXd Sigma1(3, 3);

  Sigma1 << 0.5, 0.0, 0.0,
    0.0, 0.3, 0.0,
    0.0, 0.0, 0.2;
  
  Eigen::MatrixXd Sigma2(2, 2);

  Sigma2 << 0.5, 0.0,
    0.0, 0.3;

  std::vector<Eigen::Matrix<double, -1, 1>> x1(3);
  for (size_t i = 0; i < x1.size(); ++i) {
    x1[i].resize(3, 1);
    x1[i] << 2 * i, 3 * i, 4 * i;
  }

  std::vector<Eigen::Matrix<double, -1, 1>> x2(4);
  for (size_t i = 0; i < x2.size(); ++i) {
    x2[i].resize(2, 1);
    x2[i] << 1 * i, 2 * i;
  }

  EXPECT_THROW(stan::math::gp_dot_prod_cov(x1, x2, Sigma1), std::invalid_argument);
  EXPECT_THROW(stan::math::gp_dot_prod_cov(x2, x1, Sigma1), std::invalid_argument);
  EXPECT_THROW(stan::math::gp_dot_prod_cov(x1, x1, Sigma2), std::invalid_argument);
}

TEST(MathPrimMat, vec_vec_x1_x2_gp_dot_prod_cov0) {
  Eigen::MatrixXd Sigma(3, 3);

  Sigma << 1.1, 0.30, 0.10,
    0.3, 2.70, 0.25,
    0.1, 0.25, 4.70;

  Eigen::MatrixXd ref(3, 4);

  ref << -70.45, -133.1, -195.75, -258.4,
    -140.90, -266.2, -391.50, -516.8,
    -211.35, -399.3, -587.25, -775.2;

  std::vector<Eigen::Matrix<double, -1, 1>> x1(3);
  for (size_t i = 0; i < x1.size(); ++i) {
    x1[i].resize(3, 1);
    x1[i] << 2 * (i + 1), 3 * (i + 1), 4 * (i + 1);
  }

  std::vector<Eigen::Matrix<double, -1, 1>> x2(4);
  for (size_t i = 0; i < x2.size(); ++i) {
    x2[i].resize(3, 1);
    x2[i] << 2 * (i + 1.0) - 5, 3 * (i + 1.0) + 1, -5 * (i + 1.0);
  }

  Eigen::MatrixXd cov;
  EXPECT_NO_THROW(cov = stan::math::gp_dot_prod_cov(x1, x2, Sigma));
  EXPECT_EQ(ref.rows(), cov.rows());
  EXPECT_EQ(ref.cols(), cov.cols());
  for (int i = 0; i < ref.rows(); ++i) {
    for (int j = 0; j < ref.cols(); ++j) {
      EXPECT_FLOAT_EQ(ref(i, j), cov(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }
}
