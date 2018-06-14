#include <gtest/gtest.h>
#include <limits>
#include <stan/math/prim/mat.hpp>
#include <stan/math/prim/mat/fun/squared_distance.hpp>
#include <string>
#include <vector>

template <typename T_x1, typename T_x2, typename T_sigma, typename T_l>
std::string pull_msg(std::vector<T_x1> x1, std::vector<T_x2> x2, T_sigma sigma,
                     T_l l) {
  std::string message;
  try {
    stan::math::cov_exp_quad(x1, x2, sigma, l);
  } catch (std::domain_error &e) {
    message = e.what();
  } catch (...) {
    message = "Threw the wrong exception";
  }
  return message;
}

template <typename T_x1, typename T_sigma, typename T_l>
std::string pull_msg(std::vector<T_x1> x1, T_sigma sigma, T_l l) {
  std::string message;
  try {
    stan::math::cov_exp_quad(x1, sigma, l);
  } catch (std::domain_error &e) {
    message = e.what();
  } catch (...) {
    message = "Threw the wrong exception";
  }
  return message;
}

// current proper ARD test case
TEST(MathPrimMat, ard_eigen_mat_double_cov_exp_quad1) {
  double sigma = 0.2;

  std::vector<double> l(3);
  for (int i = 0; i < 3; ++i) {
    l[i] = i + 1;
  }

  std::vector<Eigen::Matrix<double, -1, 1>> x(3);
  for (size_t i = 0; i < x.size(); ++i) {
    x[i].resize(3, 1);
    x[i] << 1 * i, 2 * i, 3 * i;
  }

  Eigen::MatrixXd cov;
  EXPECT_NO_THROW(cov = stan::math::cov_exp_quad(x, sigma, l));

  std::vector<Eigen::Matrix<double, -1, 1>> x_new(3);
  for (size_t i = 0; i < x.size(); ++i) {
    x[i].resize(3, 1);
    x[i] << 1 * i / l[0], 2 * i / l[1], 3 * i / l[2];
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_FLOAT_EQ(sigma * sigma *
                          exp(-.5 * stan::math::squared_distance(x[i], x[j])),
                      cov(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }
}
///////////

TEST(MathPrimMat, vec_double_cov_exp_quad1) {
  double sigma = 0.2;
  double l = 5;

  std::vector<double> x(3);
  x[0] = -2;
  x[1] = -1;
  x[2] = -0.5;

  Eigen::MatrixXd cov;
  cov = stan::math::cov_exp_quad(x, sigma, l);
  EXPECT_NO_THROW(cov = stan::math::cov_exp_quad(x, sigma, l));

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      EXPECT_FLOAT_EQ(sigma * sigma * exp(pow(x[i] - x[j], 2) / (-2.0 * l * l)),
                      cov(i, j))
          << "index: (" << i << ", " << j << ")";
}

TEST(MathPrimMat, vec_eigen_cov_exp_quad1) {
  using stan::math::squared_distance;
  double sigma = 0.2;
  double l = 5;

  std::vector<Eigen::Matrix<double, -1, 1>> x(3);
  for (size_t i = 0; i < x.size(); ++i) {
    x[i].resize(3, 1);
    x[i] << 1 * i, 2 * i, 3 * i;
  }

  Eigen::MatrixXd cov;
  cov = stan::math::cov_exp_quad(x, sigma, l);
  EXPECT_NO_THROW(stan::math::cov_exp_quad(x, sigma, l));
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      EXPECT_FLOAT_EQ(sigma * sigma *
                          exp(squared_distance(x[i], x[j]) / (-2.0 * l * l)),
                      cov(i, j))
          << "index: (" << i << ", " << j << ")";
}

TEST(MathPrimMat, rvec_eigen_cov_exp_quad1) {
  using stan::math::squared_distance;
  double sigma = 0.2;
  double l = 5;

  std::vector<Eigen::Matrix<double, 1, -1>> x(3);
  for (size_t i = 0; i < x.size(); ++i) {
    x[i].resize(1, 3);
    x[i] << 1 * i, 2 * i, 3 * i;
  }

  Eigen::MatrixXd cov;
  cov = stan::math::cov_exp_quad(x, sigma, l);
  EXPECT_NO_THROW(stan::math::cov_exp_quad(x, sigma, l));
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      EXPECT_FLOAT_EQ(sigma * sigma *
                          exp(squared_distance(x[i], x[j]) / (-2.0 * l * l)),
                      cov(i, j))
          << "index: (" << i << ", " << j << ")";
}

TEST(MathPrimMat, vec_double_cov_exp_quad2) {
  double sigma = 0.2;
  double l = 5;

  std::vector<double> x1(3);
  std::vector<double> x2(4);
  x1[0] = -2;
  x1[1] = -1;
  x1[2] = -0.5;

  x2[0] = 5;
  x2[1] = 0;
  x2[2] = -4;
  x2[3] = 1.1;

  Eigen::MatrixXd cov;
  EXPECT_NO_THROW(cov = stan::math::cov_exp_quad(x1, x2, sigma, l));
  EXPECT_EQ(3, cov.rows());
  EXPECT_EQ(4, cov.cols());
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 4; j++)
      EXPECT_FLOAT_EQ(sigma * sigma *
                          exp(pow(x1[i] - x2[j], 2) / (-2.0 * l * l)),
                      cov(i, j))
          << "index: (" << i << ", " << j << ")";
}

TEST(MathPrimMat, vec_eigen_rvec_cov_exp_quad2) {
  using stan::math::squared_distance;
  double sigma = 0.2;
  double l = 5;

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
  EXPECT_NO_THROW(cov = stan::math::cov_exp_quad(x1, x2, sigma, l));
  EXPECT_EQ(3, cov.rows());
  EXPECT_EQ(4, cov.cols());
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 4; j++)
      EXPECT_FLOAT_EQ(sigma * sigma *
                          exp(squared_distance(x1[i], x2[j]) / (-2.0 * l * l)),
                      cov(i, j))
          << "index: (" << i << ", " << j << ")";

  Eigen::MatrixXd cov2;
  cov2 = stan::math::cov_exp_quad(x2, x1, sigma, l);
  EXPECT_EQ(4, cov2.rows());
  EXPECT_EQ(3, cov2.cols());
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 3; j++) {
      EXPECT_FLOAT_EQ(sigma * sigma *
                          exp(squared_distance(x2[i], x1[j]) / (-2.0 * l * l)),
                      cov2(i, j))
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(cov2(i, j), cov(j, i));
    }
}

TEST(MathPrimMat, vec_eigen_vec_cov_exp_quad2) {
  using stan::math::squared_distance;
  double sigma = 0.2;
  double l = 5;

  std::vector<Eigen::Matrix<double, -1, 1>> x1(3);
  std::vector<Eigen::Matrix<double, -1, 1>> x2(4);

  for (size_t i = 0; i < x1.size(); ++i) {
    x1[i].resize(3, 1);
    x1[i] << 1 * i, 2 * i, 3 * i;
  }

  for (size_t i = 0; i < x2.size(); ++i) {
    x2[i].resize(3, 1);
    x2[i] << 2 * i, 3 * i, 4 * i;
  }

  Eigen::MatrixXd cov;
  EXPECT_NO_THROW(cov = stan::math::cov_exp_quad(x1, x2, sigma, l));
  EXPECT_EQ(3, cov.rows());
  EXPECT_EQ(4, cov.cols());
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 4; j++)
      EXPECT_FLOAT_EQ(sigma * sigma *
                          exp(squared_distance(x1[i], x2[j]) / (-2.0 * l * l)),
                      cov(i, j))
          << "index: (" << i << ", " << j << ")";

  Eigen::MatrixXd cov2;
  cov2 = stan::math::cov_exp_quad(x2, x1, sigma, l);
  EXPECT_EQ(4, cov2.rows());
  EXPECT_EQ(3, cov2.cols());
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 3; j++) {
      EXPECT_FLOAT_EQ(sigma * sigma *
                          exp(squared_distance(x2[i], x1[j]) / (-2.0 * l * l)),
                      cov2(i, j))
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(cov2(i, j), cov(j, i));
    }
}

TEST(MathPrimMat, vec_eigen_mixed_cov_exp_quad2) {
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
  EXPECT_NO_THROW(cov = stan::math::cov_exp_quad(x1_rvec, x2_vec, sigma, l));
  EXPECT_EQ(3, cov.rows());
  EXPECT_EQ(4, cov.cols());
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 4; j++)
      EXPECT_FLOAT_EQ(
          sigma * sigma *
              exp(squared_distance(x1_rvec[i], x2_vec[j]) / (-2.0 * l * l)),
          cov(i, j))
          << "index: (" << i << ", " << j << ")";

  Eigen::MatrixXd cov7;
  EXPECT_NO_THROW(cov7 = stan::math::cov_exp_quad(x2_vec, x1_rvec, sigma, l));
  EXPECT_EQ(4, cov7.rows());
  EXPECT_EQ(3, cov7.cols());
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 3; j++) {
      EXPECT_FLOAT_EQ(
          sigma * sigma *
              exp(squared_distance(x2_vec[i], x1_rvec[j]) / (-2.0 * l * l)),
          cov7(i, j))
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(cov7(i, j), cov(j, i));
    }

  Eigen::MatrixXd cov2;
  EXPECT_NO_THROW(cov2 = stan::math::cov_exp_quad(x1_vec, x2_rvec, sigma, l));
  EXPECT_EQ(3, cov2.rows());
  EXPECT_EQ(4, cov2.cols());
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 4; j++)
      EXPECT_FLOAT_EQ(
          sigma * sigma *
              exp(squared_distance(x1_vec[i], x2_rvec[j]) / (-2.0 * l * l)),
          cov2(i, j))
          << "index: (" << i << ", " << j << ")";

  Eigen::MatrixXd cov8;
  EXPECT_NO_THROW(cov8 = stan::math::cov_exp_quad(x2_rvec, x1_vec, sigma, l));
  EXPECT_EQ(4, cov8.rows());
  EXPECT_EQ(3, cov8.cols());
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 3; j++) {
      EXPECT_FLOAT_EQ(
          sigma * sigma *
              exp(squared_distance(x2_rvec[i], x1_vec[j]) / (-2.0 * l * l)),
          cov8(i, j))
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(cov8(i, j), cov2(j, i));
    }

  Eigen::MatrixXd cov3;
  EXPECT_NO_THROW(cov3 = stan::math::cov_exp_quad(x2_vec, x2_rvec, sigma, l));
  EXPECT_EQ(4, cov3.rows());
  EXPECT_EQ(4, cov3.cols());
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
      EXPECT_FLOAT_EQ(
          sigma * sigma *
              exp(squared_distance(x2_vec[i], x2_rvec[j]) / (-2.0 * l * l)),
          cov3(i, j))
          << "index: (" << i << ", " << j << ")";

  Eigen::MatrixXd cov4;
  EXPECT_NO_THROW(cov4 = stan::math::cov_exp_quad(x2_rvec, x2_vec, sigma, l));
  EXPECT_EQ(4, cov4.rows());
  EXPECT_EQ(4, cov4.cols());
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++) {
      EXPECT_FLOAT_EQ(
          sigma * sigma *
              exp(squared_distance(x2_rvec[i], x2_vec[j]) / (-2.0 * l * l)),
          cov4(i, j))
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(cov4(i, j), cov3(i, j));
    }

  Eigen::MatrixXd cov5;
  EXPECT_NO_THROW(cov5 = stan::math::cov_exp_quad(x1_rvec, x1_vec, sigma, l));
  EXPECT_EQ(3, cov5.rows());
  EXPECT_EQ(3, cov5.cols());
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      EXPECT_FLOAT_EQ(
          sigma * sigma *
              exp(squared_distance(x1_rvec[i], x1_vec[j]) / (-2.0 * l * l)),
          cov5(i, j))
          << "index: (" << i << ", " << j << ")";

  Eigen::MatrixXd cov6;
  EXPECT_NO_THROW(cov6 = stan::math::cov_exp_quad(x1_vec, x1_rvec, sigma, l));
  EXPECT_EQ(3, cov6.rows());
  EXPECT_EQ(3, cov6.cols());
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      EXPECT_FLOAT_EQ(
          sigma * sigma *
              exp(squared_distance(x1_vec[i], x1_rvec[j]) / (-2.0 * l * l)),
          cov6(i, j))
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(cov6(i, j), cov5(i, j));
    }
}

TEST(MathPrimMat, domain_error_training_sig_l) {
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

  double sigma_bad = -1;
  double l_bad = -1;

  std::string msg1, msg2, msg3;
  msg1 = pull_msg(x, sigma, l_bad);
  msg2 = pull_msg(x, sigma_bad, l);
  msg3 = pull_msg(x, sigma_bad, l_bad);
  EXPECT_TRUE(std::string::npos != msg1.find(" length-scale")) << msg1;
  EXPECT_TRUE(std::string::npos != msg2.find(" marginal variance")) << msg2;
  EXPECT_TRUE(std::string::npos != msg3.find(" marginal variance")) << msg3;

  EXPECT_THROW(stan::math::cov_exp_quad(x, sigma, l_bad), std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x, sigma_bad, l), std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x, sigma_bad, l_bad),
               std::domain_error);

  EXPECT_THROW(stan::math::cov_exp_quad(x_2, sigma, l_bad), std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_2, sigma_bad, l), std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_2, sigma_bad, l_bad),
               std::domain_error);

  EXPECT_THROW(stan::math::cov_exp_quad(x_3, sigma, l_bad), std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_3, sigma_bad, l), std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_3, sigma_bad, l_bad),
               std::domain_error);
}

TEST(MathPrimMat, nan_error_training_sig_l) {
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

  std::vector<double> x_bad(x);
  x_bad[1] = std::numeric_limits<double>::quiet_NaN();

  std::vector<Eigen::Matrix<double, -1, 1>> x_bad_2(x_2);
  x_bad_2[1](1) = std::numeric_limits<double>::quiet_NaN();

  std::vector<Eigen::Matrix<double, 1, -1>> x_bad_3(x_3);
  x_bad_3[1](1) = std::numeric_limits<double>::quiet_NaN();

  double sigma_bad = std::numeric_limits<double>::quiet_NaN();
  double l_bad = std::numeric_limits<double>::quiet_NaN();

  std::string msg1, msg2, msg3;
  msg1 = pull_msg(x, sigma, l_bad);
  msg2 = pull_msg(x, sigma_bad, l);
  msg3 = pull_msg(x, sigma_bad, l_bad);
  EXPECT_TRUE(std::string::npos != msg1.find(" length-scale")) << msg1;
  EXPECT_TRUE(std::string::npos != msg2.find(" marginal variance")) << msg2;
  EXPECT_TRUE(std::string::npos != msg3.find(" marginal variance")) << msg3;

  EXPECT_THROW(stan::math::cov_exp_quad(x, sigma, l_bad), std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x, sigma_bad, l), std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x, sigma_bad, l_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_bad, sigma, l), std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_bad, sigma_bad, l),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_bad, sigma, l_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_bad, sigma_bad, l_bad),
               std::domain_error);

  EXPECT_THROW(stan::math::cov_exp_quad(x_2, sigma, l_bad), std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_2, sigma_bad, l), std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_2, sigma_bad, l_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_bad_2, sigma, l), std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_bad_2, sigma_bad, l),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_bad_2, sigma, l_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_bad_2, sigma_bad, l_bad),
               std::domain_error);

  EXPECT_THROW(stan::math::cov_exp_quad(x_3, sigma, l_bad), std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_3, sigma_bad, l), std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_3, sigma_bad, l_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_bad_3, sigma, l), std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_bad_3, sigma_bad, l),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_bad_3, sigma, l_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_bad_3, sigma_bad, l_bad),
               std::domain_error);
}

TEST(MathPrimMat, domain_error_cov_exp_quad2) {
  double sigma = 0.2;
  double l = 5;

  std::vector<double> x1(3);
  x1[0] = -2;
  x1[1] = -1;
  x1[2] = -0.5;

  std::vector<double> x2(4);
  x2[0] = -2;
  x2[1] = -1;
  x2[2] = -0.5;
  x2[3] = -5;

  double sigma_bad = -1;
  double l_bad = -1;

  std::string msg1, msg2, msg3;
  msg1 = pull_msg(x1, x2, sigma, l_bad);
  msg2 = pull_msg(x1, x2, sigma_bad, l);
  msg3 = pull_msg(x1, x2, sigma_bad, l_bad);
  EXPECT_TRUE(std::string::npos != msg1.find(" length-scale")) << msg1;
  EXPECT_TRUE(std::string::npos != msg2.find(" marginal variance")) << msg2;
  EXPECT_TRUE(std::string::npos != msg3.find(" marginal variance")) << msg3;

  EXPECT_THROW(stan::math::cov_exp_quad(x1, x2, sigma, l_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x1, x2, sigma_bad, l),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x1, x2, sigma_bad, l_bad),
               std::domain_error);

  std::vector<Eigen::Matrix<double, -1, 1>> x_vec_1(3);
  for (size_t i = 0; i < x_vec_1.size(); ++i) {
    x_vec_1[i].resize(3, 1);
    x_vec_1[i] << 1, 2, 3;
  }

  std::vector<Eigen::Matrix<double, -1, 1>> x_vec_2(4);
  for (size_t i = 0; i < x_vec_2.size(); ++i) {
    x_vec_2[i].resize(3, 1);
    x_vec_2[i] << 4, 1, 3;
  }

  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1, x_vec_2, sigma_bad, l),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1, x_vec_2, sigma, l_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1, x_vec_2, sigma_bad, l_bad),
               std::domain_error);

  std::vector<Eigen::Matrix<double, 1, -1>> x_rvec_1(3);
  for (size_t i = 0; i < x_rvec_1.size(); ++i) {
    x_rvec_1[i].resize(1, 3);
    x_rvec_1[i] << 1, 2, 3;
  }

  std::vector<Eigen::Matrix<double, 1, -1>> x_rvec_2(4);
  for (size_t i = 0; i < x_rvec_2.size(); ++i) {
    x_rvec_2[i].resize(1, 3);
    x_rvec_2[i] << 4, 1, 3;
  }

  EXPECT_THROW(stan::math::cov_exp_quad(x_rvec_1, x_rvec_2, sigma_bad, l),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_rvec_1, x_rvec_2, sigma, l_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_rvec_1, x_rvec_2, sigma_bad, l_bad),
               std::domain_error);

  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1, x_rvec_2, sigma_bad, l),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1, x_rvec_2, sigma, l_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1, x_rvec_2, sigma_bad, l_bad),
               std::domain_error);

  EXPECT_THROW(stan::math::cov_exp_quad(x_rvec_1, x_vec_2, sigma_bad, l),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_rvec_1, x_vec_2, sigma, l_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_rvec_1, x_vec_2, sigma_bad, l_bad),
               std::domain_error);
}

TEST(MathPrimMat, nan_domain_error_cov_exp_quad2) {
  double sigma = 0.2;
  double l = 5;

  std::vector<double> x1(3);
  x1[0] = -2;
  x1[1] = -1;
  x1[2] = -0.5;

  std::vector<double> x2(4);
  x2[0] = -2;
  x2[1] = -1;
  x2[2] = -0.5;
  x2[3] = -5;

  double sigma_bad = std::numeric_limits<double>::quiet_NaN();
  double l_bad = std::numeric_limits<double>::quiet_NaN();

  std::string msg1, msg2, msg3;
  msg1 = pull_msg(x1, x2, sigma, l_bad);
  msg2 = pull_msg(x1, x2, sigma_bad, l);
  msg3 = pull_msg(x1, x2, sigma_bad, l_bad);
  EXPECT_TRUE(std::string::npos != msg1.find(" length-scale")) << msg1;
  EXPECT_TRUE(std::string::npos != msg2.find(" marginal variance")) << msg2;
  EXPECT_TRUE(std::string::npos != msg3.find(" marginal variance")) << msg3;

  EXPECT_THROW(stan::math::cov_exp_quad(x1, x2, sigma, l_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x1, x2, sigma_bad, l),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x1, x2, sigma_bad, l_bad),
               std::domain_error);

  std::vector<Eigen::Matrix<double, -1, 1>> x_vec_1(3);
  for (size_t i = 0; i < x_vec_1.size(); ++i) {
    x_vec_1[i].resize(3, 1);
    x_vec_1[i] << 1, 2, 3;
  }

  std::vector<Eigen::Matrix<double, -1, 1>> x_vec_2(4);
  for (size_t i = 0; i < x_vec_2.size(); ++i) {
    x_vec_2[i].resize(3, 1);
    x_vec_2[i] << 4, 1, 3;
  }

  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1, x_vec_2, sigma_bad, l),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1, x_vec_2, sigma, l_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1, x_vec_2, sigma_bad, l_bad),
               std::domain_error);

  std::vector<Eigen::Matrix<double, 1, -1>> x_rvec_1(3);
  for (size_t i = 0; i < x_rvec_1.size(); ++i) {
    x_rvec_1[i].resize(1, 3);
    x_rvec_1[i] << 1, 2, 3;
  }

  std::vector<Eigen::Matrix<double, 1, -1>> x_rvec_2(4);
  for (size_t i = 0; i < x_rvec_2.size(); ++i) {
    x_rvec_2[i].resize(1, 3);
    x_rvec_2[i] << 4, 1, 3;
  }

  EXPECT_THROW(stan::math::cov_exp_quad(x_rvec_1, x_rvec_2, sigma_bad, l),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_rvec_1, x_rvec_2, sigma, l_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_rvec_1, x_rvec_2, sigma_bad, l_bad),
               std::domain_error);

  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1, x_rvec_2, sigma_bad, l),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1, x_rvec_2, sigma, l_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1, x_rvec_2, sigma_bad, l_bad),
               std::domain_error);

  EXPECT_THROW(stan::math::cov_exp_quad(x_rvec_1, x_vec_2, sigma_bad, l),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_rvec_1, x_vec_2, sigma, l_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_rvec_1, x_vec_2, sigma_bad, l_bad),
               std::domain_error);

  std::vector<double> x1_bad(x1);
  x1_bad[1] = std::numeric_limits<double>::quiet_NaN();
  std::vector<double> x2_bad(x2);
  x2_bad[1] = std::numeric_limits<double>::quiet_NaN();

  std::vector<Eigen::Matrix<double, -1, 1>> x_vec_1_bad(x_vec_1);
  x_vec_1_bad[1](1) = std::numeric_limits<double>::quiet_NaN();
  std::vector<Eigen::Matrix<double, -1, 1>> x_vec_2_bad(x_vec_2);
  x_vec_2_bad[1](1) = std::numeric_limits<double>::quiet_NaN();

  std::vector<Eigen::Matrix<double, 1, -1>> x_rvec_1_bad(x_rvec_1);
  x_rvec_1_bad[1](1) = std::numeric_limits<double>::quiet_NaN();
  std::vector<Eigen::Matrix<double, 1, -1>> x_rvec_2_bad(x_rvec_2);
  x_rvec_2_bad[1](1) = std::numeric_limits<double>::quiet_NaN();

  EXPECT_THROW(stan::math::cov_exp_quad(x1_bad, x2, sigma, l),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x1, x2_bad, sigma, l),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x1_bad, x2_bad, sigma, l),
               std::domain_error);

  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1_bad, x_vec_2, sigma, l),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1, x_vec_2_bad, sigma, l),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1_bad, x_vec_2_bad, sigma, l),
               std::domain_error);

  EXPECT_THROW(stan::math::cov_exp_quad(x_rvec_1_bad, x_rvec_2, sigma, l),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_rvec_1, x_rvec_2_bad, sigma, l),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_rvec_1_bad, x_rvec_2_bad, sigma, l),
               std::domain_error);

  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1_bad, x_rvec_2, sigma, l),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_rvec_1_bad, x_vec_2, sigma, l),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1, x_rvec_2_bad, sigma, l),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_rvec_1, x_vec_2_bad, sigma, l),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1_bad, x_rvec_2_bad, sigma, l),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_rvec_1_bad, x_vec_2_bad, sigma, l),
               std::domain_error);
}

TEST(MathPrimMat, dim_mismatch_vec_eigen_vec_cov_exp_quad2) {
  double sigma = 0.2;
  double l = 5;

  std::vector<Eigen::Matrix<double, -1, 1>> x_vec_1(3);
  for (size_t i = 0; i < x_vec_1.size(); ++i) {
    x_vec_1[i].resize(3, 1);
    x_vec_1[i] << 1, 2, 3;
  }

  std::vector<Eigen::Matrix<double, -1, 1>> x_vec_2(4);
  for (size_t i = 0; i < x_vec_2.size(); ++i) {
    x_vec_2[i].resize(4, 1);
    x_vec_2[i] << 4, 1, 3, 1;
  }
  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1, x_vec_2, sigma, l),
               std::invalid_argument);
}

TEST(MathPrimMat, dim_mismatch_vec_eigen_rvec_cov_exp_quad2) {
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
  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1, x_vec_2, sigma, l),
               std::invalid_argument);
}

TEST(MathPrimMat, dim_mismatch_vec_eigen_mixed_cov_exp_quad2) {
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
  EXPECT_THROW(stan::math::cov_exp_quad(x_rvec_1, x_vec_1, sigma, l),
               std::invalid_argument);
  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1, x_rvec_1, sigma, l),
               std::invalid_argument);
  EXPECT_THROW(stan::math::cov_exp_quad(x_rvec_2, x_vec_2, sigma, l),
               std::invalid_argument);
  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_2, x_rvec_2, sigma, l),
               std::invalid_argument);
}

template <typename T_x1, typename T_x2, typename T_sigma, typename T_l>
std::string pull_msg(std::vector<T_x1> x1, std::vector<T_x2> x2, T_sigma sigma,
                     std::vector<T_l> l) {
  std::string message;
  try {
    stan::math::cov_exp_quad(x1, x2, sigma, l);
  } catch (std::domain_error &e) {
    message = e.what();
  } catch (...) {
    message = "Threw the wrong exception";
  }
  return message;
}

template <typename T_x1, typename T_sigma, typename T_l>
std::string pull_msg(std::vector<T_x1> x1, T_sigma sigma, std::vector<T_l> l) {
  std::string message;
  try {
    stan::math::cov_exp_quad(x1, sigma, l);
  } catch (std::domain_error &e) {
    message = e.what();
  } catch (...) {
    message = "Threw the wrong exception";
  }
  return message;
}

TEST(MathPrimMat, vec_length_scale_eigen_cov_exp_quad1) {
  double sigma = 0.2;

  std::vector<double> l(3);
  for (int i = 0; i < 3; ++i) {
    l[i] = i + 1;
  }

  std::vector<Eigen::Matrix<double, -1, 1>> x(3);
  for (size_t i = 0; i < x.size(); ++i) {
    x[i].resize(3, 1);
    x[i] << 1 * i, 2 * i, 3 * i;
  }

  Eigen::MatrixXd cov;
  EXPECT_NO_THROW(cov = stan::math::cov_exp_quad(x, sigma, l));

  std::vector<Eigen::Matrix<double, -1, 1>> x_new(3);
  for (size_t i = 0; i < x.size(); ++i) {
    x_new[i].resize(3, 1);
    x_new[i] << 1 * i / l[0], 2 * i / l[1], 3 * i / l[2];
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_FLOAT_EQ(sigma * sigma * exp(-.5 * stan::math::squared_distance(
                                                    x_new[i], x_new[j])),
                      cov(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }
}

TEST(MathPrimMat, rvec_length_scale_eigen_cov_exp_quad1) {
  using stan::math::squared_distance;
  double sigma = 0.2;

  std::vector<Eigen::Matrix<double, 1, -1>> x(3);
  for (size_t i = 0; i < x.size(); ++i) {
    x[i].resize(1, 3);
    x[i] << 1 * i, 2 * i, 3 * i;
  }

  std::vector<double> l(3);
  for (int i = 0; i < 3; ++i) {
    l[i] = i + 1;
  }

  Eigen::MatrixXd cov;
  EXPECT_NO_THROW(cov = stan::math::cov_exp_quad(x, sigma, l));
  std::vector<Eigen::Matrix<double, 1, -1>> x_new(3);
  for (size_t i = 0; i < 3; ++i) {
    x_new[i].resize(1, 3);
    x_new[i] << 1 * i / l[0], 2 * i / l[1], 3 * i / l[2];
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_FLOAT_EQ(sigma * sigma * exp(-.5 * stan::math::squared_distance(
                                                    x_new[i], x_new[j])),
                      cov(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }
}

TEST(MathPrimMat, vec_length_scale_eigen_rvec_cov_exp_quad2) {
  using stan::math::squared_distance;
  double sigma = 0.2;

  std::vector<double> l(3);
  for (int i = 0; i < 3; ++i) {
    l[i] = i + 1;
  }

  std::vector<Eigen::Matrix<double, 1, -1>> x1(3);
  std::vector<Eigen::Matrix<double, 1, -1>> x2(3);

  for (size_t i = 0; i < x1.size(); ++i) {
    x1[i].resize(1, 3);
    x1[i] << 1 * i, 2 * i, 3 * i;
  }

  for (size_t i = 0; i < x2.size(); ++i) {
    x2[i].resize(1, 3);
    x2[i] << 2 * i, 3 * i, 4 * i;
  }

  std::vector<Eigen::Matrix<double, 1, -1>> x1_new(3);
  for (size_t i = 0; i < 3; ++i) {
    x1_new[i].resize(1, 3);
    x1_new[i] << 1 * i / l[0], 2 * i / l[1], 3 * i / l[2];
  }

  std::vector<Eigen::Matrix<double, 1, -1>> x2_new(3);
  for (size_t i = 0; i < 3; ++i) {
    x2_new[i].resize(1, 3);
    x2_new[i] << 2 * i / l[0], 3 * i / l[1], 4 * i / l[2];
  }

  Eigen::MatrixXd cov;
  EXPECT_NO_THROW(cov = stan::math::cov_exp_quad(x1, x2, sigma, l));
  EXPECT_EQ(3, cov.rows());
  EXPECT_EQ(3, cov.cols());
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_FLOAT_EQ(sigma * sigma * exp(-.5 * stan::math::squared_distance(
                                                    x1_new[i], x2_new[j])),
                      cov(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }

  Eigen::MatrixXd cov2;
  cov2 = stan::math::cov_exp_quad(x2, x1, sigma, l);
  EXPECT_EQ(3, cov2.rows());
  EXPECT_EQ(3, cov2.cols());
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_FLOAT_EQ(sigma * sigma * exp(-.5 * stan::math::squared_distance(
                                                    x2_new[i], x1_new[j])),
                      cov2(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }
}

TEST(MathPrimMat, vec_length_scale_vec_eigen_vec_cov_exp_quad2) {
  using stan::math::squared_distance;
  double sigma = 0.2;

  std::vector<double> l(3);
  for (int i = 0; i < 3; ++i) {
    l[i] = i + 1;
  }

  std::vector<Eigen::Matrix<double, -1, 1>> x1(3);
  std::vector<Eigen::Matrix<double, -1, 1>> x2(3);

  for (size_t i = 0; i < x1.size(); ++i) {
    x1[i].resize(3, 1);
    x1[i] << 1 * i, 2 * i, 3 * i;
  }

  for (size_t i = 0; i < x2.size(); ++i) {
    x2[i].resize(3, 1);
    x2[i] << 2 * i, 3 * i, 4 * i;
  }

  std::vector<Eigen::Matrix<double, 1, -1>> x1_new(3);
  for (size_t i = 0; i < 3; ++i) {
    x1_new[i].resize(1, 3);
    x1_new[i] << 1 * i / l[0], 2 * i / l[1], 3 * i / l[2];
  }

  std::vector<Eigen::Matrix<double, 1, -1>> x2_new(3);
  for (size_t i = 0; i < 3; ++i) {
    x2_new[i].resize(1, 3);
    x2_new[i] << 2 * i / l[0], 3 * i / l[1], 4 * i / l[2];
  }

  Eigen::MatrixXd cov;
  EXPECT_NO_THROW(cov = stan::math::cov_exp_quad(x1, x2, sigma, l));
  EXPECT_EQ(3, cov.rows());
  EXPECT_EQ(3, cov.cols());
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_FLOAT_EQ(sigma * sigma * exp(-.5 * stan::math::squared_distance(
                                                    x1_new[i], x2_new[j])),
                      cov(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }
  Eigen::MatrixXd cov2;
  cov2 = stan::math::cov_exp_quad(x2, x1, sigma, l);
  EXPECT_EQ(3, cov2.rows());
  EXPECT_EQ(3, cov2.cols());
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_FLOAT_EQ(sigma * sigma * exp(-.5 * stan::math::squared_distance(
                                                    x2_new[i], x1_new[j])),
                      cov2(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }
}

TEST(MathPrimMat, vec_length_scale_vec_eigen_mixed_cov_exp_quad2) {
  using stan::math::squared_distance;
  double sigma = 0.2;

  std::vector<double> l(3);
  for (int i = 0; i < 3; ++i) {
    l[i] = i + 1;
  }

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

  std::vector<Eigen::Matrix<double, 1, -1>> x1_rvec_new(3);
  for (size_t i = 0; i < 3; ++i) {
    x1_rvec_new[i].resize(1, 3);
    x1_rvec_new[i] << 1 * i / l[0], 2 * i / l[1], 3 * i / l[2];
  }

  std::vector<Eigen::Matrix<double, 1, -1>> x2_rvec_new(4);
  for (size_t i = 0; i < 4; ++i) {
    x2_rvec_new[i].resize(1, 3);
    x2_rvec_new[i] << 2 * i / l[0], 3 * i / l[1], 4 * i / l[2];
  }

  std::vector<Eigen::Matrix<double, -1, 1>> x1_vec_new(3);
  for (size_t i = 0; i < 3; ++i) {
    x1_vec_new[i].resize(3, 1);
    x1_vec_new[i] << 1 * i / l[0], 2 * i / l[1], 3 * i / l[2];
  }

  std::vector<Eigen::Matrix<double, -1, 1>> x2_vec_new(4);
  for (size_t i = 0; i < 4; ++i) {
    x2_vec_new[i].resize(3, 1);
    x2_vec_new[i] << 2 * i / l[0], 3 * i / l[1], 4 * i / l[2];
  }

  Eigen::MatrixXd cov;
  EXPECT_NO_THROW(cov = stan::math::cov_exp_quad(x1_rvec, x2_vec, sigma, l));
  EXPECT_EQ(3, cov.rows());
  EXPECT_EQ(4, cov.cols());
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 4; j++) {
      EXPECT_FLOAT_EQ(sigma * sigma *
                          exp(-.5 * stan::math::squared_distance(
                                        x1_rvec_new[i], x2_vec_new[j])),
                      cov(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }

  Eigen::MatrixXd cov7;
  EXPECT_NO_THROW(cov7 = stan::math::cov_exp_quad(x2_vec, x1_rvec, sigma, l));
  EXPECT_EQ(4, cov7.rows());
  EXPECT_EQ(3, cov7.cols());
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_FLOAT_EQ(sigma * sigma *
                          exp(-.5 * stan::math::squared_distance(
                                        x2_vec_new[i], x1_rvec_new[j])),
                      cov7(i, j))
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(cov7(i, j), cov(j, i));
    }
  }

  Eigen::MatrixXd cov2;
  EXPECT_NO_THROW(cov2 = stan::math::cov_exp_quad(x1_vec, x2_rvec, sigma, l));
  EXPECT_EQ(3, cov2.rows());
  EXPECT_EQ(4, cov2.cols());
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 4; j++) {
      EXPECT_FLOAT_EQ(sigma * sigma *
                          exp(-.5 * stan::math::squared_distance(
                                        x1_vec_new[i], x2_rvec_new[j])),
                      cov2(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }

  Eigen::MatrixXd cov8;
  EXPECT_NO_THROW(cov8 = stan::math::cov_exp_quad(x2_rvec, x1_vec, sigma, l));
  EXPECT_EQ(4, cov8.rows());
  EXPECT_EQ(3, cov8.cols());
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_FLOAT_EQ(sigma * sigma *
                          exp(-.5 * stan::math::squared_distance(
                                        x2_rvec_new[i], x1_vec_new[j])),
                      cov8(i, j))
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(cov8(i, j), cov2(j, i));
    }
  }

  Eigen::MatrixXd cov3;
  EXPECT_NO_THROW(cov3 = stan::math::cov_exp_quad(x2_vec, x2_rvec, sigma, l));
  EXPECT_EQ(4, cov3.rows());
  EXPECT_EQ(4, cov3.cols());
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      EXPECT_FLOAT_EQ(sigma * sigma *
                          exp(-.5 * stan::math::squared_distance(
                                        x2_vec_new[i], x2_rvec_new[j])),
                      cov3(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }

  Eigen::MatrixXd cov4;
  EXPECT_NO_THROW(cov4 = stan::math::cov_exp_quad(x2_rvec, x2_vec, sigma, l));
  EXPECT_EQ(4, cov4.rows());
  EXPECT_EQ(4, cov4.cols());
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      EXPECT_FLOAT_EQ(sigma * sigma *
                          exp(-.5 * stan::math::squared_distance(
                                        x2_rvec_new[i], x2_vec_new[j])),
                      cov4(i, j))
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(cov4(i, j), cov3(i, j));
    }
  }

  Eigen::MatrixXd cov5;
  EXPECT_NO_THROW(cov5 = stan::math::cov_exp_quad(x1_rvec, x1_vec, sigma, l));
  EXPECT_EQ(3, cov5.rows());
  EXPECT_EQ(3, cov5.cols());
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_FLOAT_EQ(sigma * sigma *
                          exp(-.5 * stan::math::squared_distance(
                                        x1_rvec_new[i], x1_vec_new[j])),
                      cov5(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }

  Eigen::MatrixXd cov6;
  EXPECT_NO_THROW(cov6 = stan::math::cov_exp_quad(x1_vec, x1_rvec, sigma, l));
  EXPECT_EQ(3, cov6.rows());
  EXPECT_EQ(3, cov6.cols());
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_FLOAT_EQ(sigma * sigma *
                          exp(-.5 * stan::math::squared_distance(
                                        x1_vec_new[i], x1_rvec_new[j])),
                      cov6(i, j))
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(cov6(i, j), cov5(j, i));
    }
  }
}

TEST(MathPrimMat, domain_error_training_sig_vec_length_scale) {
  double sigma = 0.2;

  std::vector<double> l(3);
  for (int i = 0; i < 3; ++i) {
    l[i] = i + 1;
  }

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

  double sigma_bad = -1;
  std::vector<double> l_bad(3);
  for (int i = 0; i < 3; ++i) {
    l[i] = i + 1;
  }
  l[2] = -1;

  std::string msg1, msg2, msg3;
  msg1 = pull_msg(x_2, sigma, l_bad);
  msg2 = pull_msg(x_2, sigma_bad, l);
  msg3 = pull_msg(x_2, sigma_bad, l_bad);

  EXPECT_TRUE(std::string::npos != msg1.find(" length scale")) << msg1;
  EXPECT_TRUE(std::string::npos != msg2.find(" marginal variance")) << msg2;
  EXPECT_TRUE(std::string::npos != msg3.find(" marginal variance")) << msg3;

  EXPECT_THROW(stan::math::cov_exp_quad(x_2, sigma, l_bad), std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_2, sigma_bad, l), std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_2, sigma_bad, l_bad),
               std::domain_error);

  EXPECT_THROW(stan::math::cov_exp_quad(x_3, sigma, l_bad), std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_3, sigma_bad, l), std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_3, sigma_bad, l_bad),
               std::domain_error);
}

TEST(MathPrimMat, nan_error_training_sig_vec_length_scale) {
  double sigma = 0.2;

  std::vector<double> l(3);
  for (int i = 0; i < 3; ++i) {
    l[i] = i + 1;
  }

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

  std::vector<Eigen::Matrix<double, -1, 1>> x_bad_2(x_2);
  x_bad_2[1](1) = std::numeric_limits<double>::quiet_NaN();

  std::vector<Eigen::Matrix<double, 1, -1>> x_bad_3(x_3);
  x_bad_3[1](1) = std::numeric_limits<double>::quiet_NaN();

  double sigma_bad = std::numeric_limits<double>::quiet_NaN();

  std::vector<double> l_bad(3);
  for (int i = 0; i < 3; ++i) {
    l_bad[i] = i + 1;
  }
  l_bad[1] = std::numeric_limits<double>::quiet_NaN();

  std::string msg1, msg2, msg3;
  msg1 = pull_msg(x_2, sigma, l_bad);
  EXPECT_TRUE(std::string::npos != msg1.find(" length scale")) << msg1;

  EXPECT_THROW(stan::math::cov_exp_quad(x_2, sigma, l_bad), std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_2, sigma_bad, l_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_bad_2, sigma, l_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_bad_2, sigma_bad, l_bad),
               std::domain_error);

  EXPECT_THROW(stan::math::cov_exp_quad(x_3, sigma, l_bad), std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_3, sigma_bad, l_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_bad_3, sigma, l_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_bad_3, sigma_bad, l_bad),
               std::domain_error);
}

TEST(MathPrimMat, domain_error_cov_exp_quad2_vec_length_scale) {
  double sigma = 0.2;

  std::vector<double> l(3);
  for (int i = 0; i < 3; ++i) {
    l[i] = i + 1;
  }

  double sigma_bad = -1;
  std::vector<double> l_bad0(3);
  std::vector<double> l_bad1(3);
  std::vector<double> l_bad2(3);
  for (int i = 0; i < 3; ++i) {
    l_bad0[i] = i + 1;
    l_bad1[i] = -1;
    l_bad2[i] = i + 1;
  }
  l_bad0[0] = -10;
  l_bad2[2] = -7;

  std::vector<Eigen::Matrix<double, -1, 1>> x_vec_1(3);
  for (size_t i = 0; i < x_vec_1.size(); ++i) {
    x_vec_1[i].resize(3, 1);
    x_vec_1[i] << 1, 2, 3;
  }

  std::vector<Eigen::Matrix<double, -1, 1>> x_vec_2(4);
  for (size_t i = 0; i < x_vec_2.size(); ++i) {
    x_vec_2[i].resize(3, 1);
    x_vec_2[i] << 4, 1, 3;
  }

  std::string msg1, msg2, msg3, msg4;
  msg1 = pull_msg(x_vec_1, x_vec_2, sigma, l_bad0);
  msg2 = pull_msg(x_vec_1, x_vec_2, sigma_bad, l_bad1);
  msg3 = pull_msg(x_vec_1, x_vec_2, sigma, l_bad1);
  msg4 = pull_msg(x_vec_1, x_vec_2, sigma, l_bad2);

  EXPECT_TRUE(std::string::npos != msg1.find(" length scale")) << msg1;
  EXPECT_TRUE(std::string::npos != msg2.find(" marginal variance")) << msg2;
  EXPECT_TRUE(std::string::npos != msg3.find(" length scale")) << msg3;
  EXPECT_TRUE(std::string::npos != msg4.find(" length scale")) << msg4;

  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1, x_vec_2, sigma, l_bad2),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1, x_vec_2, sigma, l_bad0),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1, x_vec_2, sigma, l_bad1),
               std::domain_error);

  std::vector<Eigen::Matrix<double, 1, -1>> x_rvec_1(3);
  for (size_t i = 0; i < x_rvec_1.size(); ++i) {
    x_rvec_1[i].resize(1, 3);
    x_rvec_1[i] << 1, 2, 3;
  }

  std::vector<Eigen::Matrix<double, 1, -1>> x_rvec_2(4);
  for (size_t i = 0; i < x_rvec_2.size(); ++i) {
    x_rvec_2[i].resize(1, 3);
    x_rvec_2[i] << 4, 1, 3;
  }

  EXPECT_THROW(stan::math::cov_exp_quad(x_rvec_1, x_rvec_2, sigma, l_bad0),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_rvec_1, x_rvec_2, sigma, l_bad1),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_rvec_1, x_rvec_2, sigma, l_bad2),
               std::domain_error);

  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1, x_rvec_2, sigma, l_bad0),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1, x_rvec_2, sigma, l_bad1),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1, x_rvec_2, sigma, l_bad2),
               std::domain_error);

  EXPECT_THROW(stan::math::cov_exp_quad(x_rvec_1, x_vec_2, sigma, l_bad0),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_rvec_1, x_vec_2, sigma, l_bad1),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_rvec_1, x_vec_2, sigma, l_bad2),
               std::domain_error);
}

TEST(MathPrimMat, nan_domain_error_cov_exp_quad2_vec_length_scale) {
  double sigma = 0.2;

  std::vector<double> l(3);
  for (int i = 0; i < 3; ++i) {
    l[i] = i + 1;
  }

  double sigma_bad = std::numeric_limits<double>::quiet_NaN();
  std::vector<double> l_bad0(3);
  std::vector<double> l_bad1(3);
  std::vector<double> l_bad2(3);
  std::vector<double> l_bad3(3);
  std::vector<double> l_bad4(3);
  for (int i = 0; i < 3; ++i) {
    l_bad0[i] = std::numeric_limits<double>::quiet_NaN();
    l_bad1[i] = i + 1;
    l_bad2[i] = i + 1;
    l_bad3[i] = i + 1;
    l_bad4[i] = i + 1;
  }
  l_bad1[0] = std::numeric_limits<double>::quiet_NaN();
  l_bad2[1] = std::numeric_limits<double>::quiet_NaN();
  l_bad3[2] = std::numeric_limits<double>::quiet_NaN();
  l_bad4[0] = std::numeric_limits<double>::quiet_NaN();
  l_bad4[1] = std::numeric_limits<double>::quiet_NaN();

  std::vector<Eigen::Matrix<double, -1, 1>> x_vec_1(3);
  for (size_t i = 0; i < x_vec_1.size(); ++i) {
    x_vec_1[i].resize(3, 1);
    x_vec_1[i] << 1, 2, 3;
  }

  std::vector<Eigen::Matrix<double, -1, 1>> x_vec_2(4);
  for (size_t i = 0; i < x_vec_2.size(); ++i) {
    x_vec_2[i].resize(3, 1);
    x_vec_2[i] << 4, 1, 3;
  }

  std::string msg1, msg2, msg3, msg4, msg5, msg6, msg7;
  msg1 = pull_msg(x_vec_1, x_vec_2, sigma, l_bad0);
  msg2 = pull_msg(x_vec_1, x_vec_2, sigma_bad, l_bad1);
  msg3 = pull_msg(x_vec_1, x_vec_2, sigma_bad, l);
  msg4 = pull_msg(x_vec_1, x_vec_2, sigma, l_bad1);
  msg5 = pull_msg(x_vec_1, x_vec_2, sigma, l_bad2);
  msg6 = pull_msg(x_vec_1, x_vec_2, sigma, l_bad3);
  msg7 = pull_msg(x_vec_1, x_vec_2, sigma, l_bad4);

  EXPECT_TRUE(std::string::npos != msg1.find(" length scale")) << msg1;
  EXPECT_TRUE(std::string::npos != msg2.find(" marginal variance")) << msg2;
  EXPECT_TRUE(std::string::npos != msg3.find(" marginal variance")) << msg3;
  EXPECT_TRUE(std::string::npos != msg4.find(" length scale")) << msg4;
  EXPECT_TRUE(std::string::npos != msg5.find(" length scale")) << msg5;
  EXPECT_TRUE(std::string::npos != msg6.find(" length scale")) << msg6;
  EXPECT_TRUE(std::string::npos != msg7.find(" length scale")) << msg7;

  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1, x_vec_2, sigma_bad, l),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1, x_vec_2, sigma, l_bad0),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1, x_vec_2, sigma_bad, l_bad0),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1, x_vec_2, sigma, l_bad1),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1, x_vec_2, sigma, l_bad2),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1, x_vec_2, sigma, l_bad3),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1, x_vec_2, sigma, l_bad4),
               std::domain_error);

  std::vector<Eigen::Matrix<double, 1, -1>> x_rvec_1(3);
  for (size_t i = 0; i < x_rvec_1.size(); ++i) {
    x_rvec_1[i].resize(1, 3);
    x_rvec_1[i] << 1, 2, 3;
  }

  std::vector<Eigen::Matrix<double, 1, -1>> x_rvec_2(4);
  for (size_t i = 0; i < x_rvec_2.size(); ++i) {
    x_rvec_2[i].resize(1, 3);
    x_rvec_2[i] << 4, 1, 3;
  }

  EXPECT_THROW(stan::math::cov_exp_quad(x_rvec_1, x_rvec_2, sigma_bad, l),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_rvec_1, x_rvec_2, sigma, l_bad0),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_rvec_1, x_rvec_2, sigma_bad, l_bad0),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_rvec_1, x_rvec_2, sigma, l_bad1),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_rvec_1, x_rvec_2, sigma, l_bad2),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_rvec_1, x_rvec_2, sigma, l_bad3),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_rvec_1, x_rvec_2, sigma, l_bad4),
               std::domain_error);

  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1, x_rvec_2, sigma_bad, l),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1, x_rvec_2, sigma, l_bad0),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1, x_rvec_2, sigma_bad, l_bad0),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1, x_rvec_2, sigma, l_bad1),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1, x_rvec_2, sigma, l_bad2),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1, x_rvec_2, sigma, l_bad3),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1, x_rvec_2, sigma, l_bad4),
               std::domain_error);

  EXPECT_THROW(stan::math::cov_exp_quad(x_rvec_1, x_vec_2, sigma_bad, l),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_rvec_1, x_vec_2, sigma, l_bad0),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_rvec_1, x_vec_2, sigma_bad, l_bad0),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_rvec_1, x_vec_2, sigma, l_bad1),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_rvec_1, x_vec_2, sigma, l_bad2),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_rvec_1, x_vec_2, sigma, l_bad3),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_rvec_1, x_vec_2, sigma, l_bad4),
               std::domain_error);

  std::vector<Eigen::Matrix<double, -1, 1>> x_vec_1_bad(x_vec_1);
  x_vec_1_bad[1](1) = std::numeric_limits<double>::quiet_NaN();
  std::vector<Eigen::Matrix<double, -1, 1>> x_vec_2_bad(x_vec_2);
  x_vec_2_bad[1](1) = std::numeric_limits<double>::quiet_NaN();

  std::vector<Eigen::Matrix<double, 1, -1>> x_rvec_1_bad(x_rvec_1);
  x_rvec_1_bad[1](1) = std::numeric_limits<double>::quiet_NaN();
  std::vector<Eigen::Matrix<double, 1, -1>> x_rvec_2_bad(x_rvec_2);
  x_rvec_2_bad[1](1) = std::numeric_limits<double>::quiet_NaN();

  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1_bad, x_vec_2, sigma, l),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1, x_vec_2_bad, sigma, l),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1_bad, x_vec_2_bad, sigma, l),
               std::domain_error);

  EXPECT_THROW(stan::math::cov_exp_quad(x_rvec_1_bad, x_rvec_2, sigma, l),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_rvec_1, x_rvec_2_bad, sigma, l),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_rvec_1_bad, x_rvec_2_bad, sigma, l),
               std::domain_error);

  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1_bad, x_rvec_2, sigma, l),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_rvec_1_bad, x_vec_2, sigma, l),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1, x_rvec_2_bad, sigma, l),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_rvec_1, x_vec_2_bad, sigma, l),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_vec_1_bad, x_rvec_2_bad, sigma, l),
               std::domain_error);
  EXPECT_THROW(stan::math::cov_exp_quad(x_rvec_1_bad, x_vec_2_bad, sigma, l),
               std::domain_error);
}
