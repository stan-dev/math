
#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <string>
#include <vector>

TEST(MathPrimMat, ard_eigen_mat_double_gp_exp_quad_cov1) {
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
  EXPECT_NO_THROW(cov = stan::math::gp_exp_quad_cov(x, sigma, l));

  std::vector<Eigen::Matrix<double, -1, 1>> x_new(3);
  for (size_t i = 0; i < x.size(); ++i) {
    x[i].resize(3, 1);
    x[i] << 1 * i / l[0], 2 * i / l[1], 3 * i / l[2];
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_FLOAT_EQ(
          sigma * sigma * exp(-.5 * stan::math::squared_distance(x[i], x[j])),
          cov(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }
}

TEST(MathPrimMat, vec_double_gp_exp_quad_cov1) {
  double sigma = 0.2;
  double l = 5;

  std::vector<double> x(3);
  x[0] = -2;
  x[1] = -1;
  x[2] = -0.5;

  Eigen::MatrixXd cov;
  cov = stan::math::gp_exp_quad_cov(x, sigma, l);
  EXPECT_NO_THROW(cov = stan::math::gp_exp_quad_cov(x, sigma, l));

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      EXPECT_FLOAT_EQ(sigma * sigma * exp(pow(x[i] - x[j], 2) / (-2.0 * l * l)),
                      cov(i, j))
          << "index: (" << i << ", " << j << ")";
}

TEST(MathPrimMat, vec_eigen_gp_exp_quad_cov1) {
  using stan::math::squared_distance;
  using std::exp;
  double sigma = 0.2;
  double l = 5;

  std::vector<Eigen::Matrix<double, -1, 1>> x(3);
  for (size_t i = 0; i < x.size(); ++i) {
    x[i].resize(3, 1);
    x[i] << 1 * i, 2 * i, 3 * i;
  }

  Eigen::MatrixXd cov;
  cov = stan::math::gp_exp_quad_cov(x, sigma, l);
  EXPECT_NO_THROW(stan::math::gp_exp_quad_cov(x, sigma, l));
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      EXPECT_FLOAT_EQ(
          sigma * sigma * exp(squared_distance(x[i], x[j]) / (-2.0 * l * l)),
          cov(i, j))
          << "index: (" << i << ", " << j << ")";
}

TEST(MathPrimMat, rvec_eigen_gp_exp_quad_cov1) {
  using stan::math::squared_distance;
  double sigma = 0.2;
  double l = 5;

  std::vector<Eigen::Matrix<double, 1, -1>> x(3);
  for (size_t i = 0; i < x.size(); ++i) {
    x[i].resize(1, 3);
    x[i] << 1 * i, 2 * i, 3 * i;
  }

  Eigen::MatrixXd cov;
  cov = stan::math::gp_exp_quad_cov(x, sigma, l);
  EXPECT_NO_THROW(stan::math::gp_exp_quad_cov(x, sigma, l));
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      EXPECT_FLOAT_EQ(
          sigma * sigma * exp(squared_distance(x[i], x[j]) / (-2.0 * l * l)),
          cov(i, j))
          << "index: (" << i << ", " << j << ")";
}

TEST(MathPrimMat, vec_double_gp_exp_quad_cov2) {
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
  EXPECT_NO_THROW(cov = stan::math::gp_exp_quad_cov(x1, x2, sigma, l));
  EXPECT_EQ(3, cov.rows());
  EXPECT_EQ(4, cov.cols());
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 4; j++)
      EXPECT_FLOAT_EQ(
          sigma * sigma * exp(pow(x1[i] - x2[j], 2) / (-2.0 * l * l)),
          cov(i, j))
          << "index: (" << i << ", " << j << ")";
}

TEST(MathPrimMat, vec_eigen_vec_gp_exp_quad_cov2) {
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
  EXPECT_NO_THROW(cov = stan::math::gp_exp_quad_cov(x1, x2, sigma, l));
  EXPECT_EQ(3, cov.rows());
  EXPECT_EQ(4, cov.cols());
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 4; j++)
      EXPECT_FLOAT_EQ(
          sigma * sigma * exp(squared_distance(x1[i], x2[j]) / (-2.0 * l * l)),
          cov(i, j))
          << "index: (" << i << ", " << j << ")";

  Eigen::MatrixXd cov2;
  cov2 = stan::math::gp_exp_quad_cov(x2, x1, sigma, l);
  EXPECT_EQ(4, cov2.rows());
  EXPECT_EQ(3, cov2.cols());
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 3; j++) {
      EXPECT_FLOAT_EQ(
          sigma * sigma * exp(squared_distance(x2[i], x1[j]) / (-2.0 * l * l)),
          cov2(i, j))
          << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(cov2(i, j), cov(j, i));
    }
}

TEST(MathPrimMat, domain_error_training_sig_l_gp_cov_exp) {
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

  EXPECT_THROW_MSG(stan::math::gp_exp_quad_cov(x, sigma, l_bad),
                   std::domain_error, " length scale");
  EXPECT_THROW_MSG(stan::math::gp_exp_quad_cov(x, sigma_bad, l),
                   std::domain_error, " magnitude");

  EXPECT_THROW(stan::math::gp_exp_quad_cov(x_2, sigma, l_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_exp_quad_cov(x_2, sigma_bad, l),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_exp_quad_cov(x_2, sigma_bad, l_bad),
               std::domain_error);

  EXPECT_THROW(stan::math::gp_exp_quad_cov(x_3, sigma, l_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_exp_quad_cov(x_3, sigma_bad, l),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_exp_quad_cov(x_3, sigma_bad, l_bad),
               std::domain_error);
}

TEST(MathPrimMat, nan_error_training_sig_l_gp_cov_exp) {
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

  EXPECT_THROW_MSG(stan::math::gp_exp_quad_cov(x, sigma, l_bad),
                   std::domain_error, " length scale");
  EXPECT_THROW_MSG(stan::math::gp_exp_quad_cov(x, sigma_bad, l),
                   std::domain_error, " magnitude");

  EXPECT_THROW(stan::math::gp_exp_quad_cov(x, sigma, l_bad), std::domain_error);
  EXPECT_THROW(stan::math::gp_exp_quad_cov(x, sigma_bad, l), std::domain_error);
  EXPECT_THROW(stan::math::gp_exp_quad_cov(x, sigma_bad, l_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_exp_quad_cov(x_bad, sigma, l), std::domain_error);
  EXPECT_THROW(stan::math::gp_exp_quad_cov(x_bad, sigma_bad, l),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_exp_quad_cov(x_bad, sigma, l_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_exp_quad_cov(x_bad, sigma_bad, l_bad),
               std::domain_error);

  EXPECT_THROW(stan::math::gp_exp_quad_cov(x_2, sigma, l_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_exp_quad_cov(x_2, sigma_bad, l),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_exp_quad_cov(x_2, sigma_bad, l_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_exp_quad_cov(x_bad_2, sigma, l),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_exp_quad_cov(x_bad_2, sigma_bad, l),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_exp_quad_cov(x_bad_2, sigma, l_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_exp_quad_cov(x_bad_2, sigma_bad, l_bad),
               std::domain_error);

  EXPECT_THROW(stan::math::gp_exp_quad_cov(x_3, sigma, l_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_exp_quad_cov(x_3, sigma_bad, l),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_exp_quad_cov(x_3, sigma_bad, l_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_exp_quad_cov(x_bad_3, sigma, l),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_exp_quad_cov(x_bad_3, sigma_bad, l),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_exp_quad_cov(x_bad_3, sigma, l_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_exp_quad_cov(x_bad_3, sigma_bad, l_bad),
               std::domain_error);
}

TEST(MathPrimMat, domain_error_gp_exp_quad_cov2) {
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

  EXPECT_THROW_MSG(stan::math::gp_exp_quad_cov(x1, x2, sigma, l_bad),
                   std::domain_error, " length scale");
  EXPECT_THROW_MSG(stan::math::gp_exp_quad_cov(x1, x2, sigma_bad, l),
                   std::domain_error, " magnitude");
  EXPECT_THROW(stan::math::gp_exp_quad_cov(x1, x2, sigma_bad, l_bad),
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

  EXPECT_THROW(stan::math::gp_exp_quad_cov(x_vec_1, x_vec_2, sigma_bad, l),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_exp_quad_cov(x_vec_1, x_vec_2, sigma, l_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_exp_quad_cov(x_vec_1, x_vec_2, sigma_bad, l_bad),
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

  EXPECT_THROW(stan::math::gp_exp_quad_cov(x_rvec_1, x_rvec_2, sigma_bad, l),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_exp_quad_cov(x_rvec_1, x_rvec_2, sigma, l_bad),
               std::domain_error);
  EXPECT_THROW(
      stan::math::gp_exp_quad_cov(x_rvec_1, x_rvec_2, sigma_bad, l_bad),
      std::domain_error);

  EXPECT_THROW(stan::math::gp_exp_quad_cov(x_vec_1, x_rvec_2, sigma_bad, l),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_exp_quad_cov(x_vec_1, x_rvec_2, sigma, l_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_exp_quad_cov(x_vec_1, x_rvec_2, sigma_bad, l_bad),
               std::domain_error);

  EXPECT_THROW(stan::math::gp_exp_quad_cov(x_rvec_1, x_vec_2, sigma_bad, l),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_exp_quad_cov(x_rvec_1, x_vec_2, sigma, l_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_exp_quad_cov(x_rvec_1, x_vec_2, sigma_bad, l_bad),
               std::domain_error);
}

TEST(MathPrimMat, nan_domain_error_gp_exp_quad_cov2) {
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

  EXPECT_THROW_MSG(stan::math::gp_exp_quad_cov(x1, x2, sigma, l_bad),
                   std::domain_error, " length scale");
  EXPECT_THROW_MSG(stan::math::gp_exp_quad_cov(x1, x2, sigma_bad, l),
                   std::domain_error, " magnitude");
  EXPECT_THROW(stan::math::gp_exp_quad_cov(x1, x2, sigma_bad, l_bad),
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

  EXPECT_THROW(stan::math::gp_exp_quad_cov(x_vec_1, x_vec_2, sigma_bad, l),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_exp_quad_cov(x_vec_1, x_vec_2, sigma, l_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_exp_quad_cov(x_vec_1, x_vec_2, sigma_bad, l_bad),
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

  EXPECT_THROW(stan::math::gp_exp_quad_cov(x_rvec_1, x_rvec_2, sigma_bad, l),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_exp_quad_cov(x_rvec_1, x_rvec_2, sigma, l_bad),
               std::domain_error);
  EXPECT_THROW(
      stan::math::gp_exp_quad_cov(x_rvec_1, x_rvec_2, sigma_bad, l_bad),
      std::domain_error);

  EXPECT_THROW(stan::math::gp_exp_quad_cov(x_vec_1, x_rvec_2, sigma_bad, l),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_exp_quad_cov(x_vec_1, x_rvec_2, sigma, l_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_exp_quad_cov(x_vec_1, x_rvec_2, sigma_bad, l_bad),
               std::domain_error);

  EXPECT_THROW(stan::math::gp_exp_quad_cov(x_rvec_1, x_vec_2, sigma_bad, l),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_exp_quad_cov(x_rvec_1, x_vec_2, sigma, l_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_exp_quad_cov(x_rvec_1, x_vec_2, sigma_bad, l_bad),
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

  EXPECT_THROW(stan::math::gp_exp_quad_cov(x1_bad, x2, sigma, l),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_exp_quad_cov(x1, x2_bad, sigma, l),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_exp_quad_cov(x1_bad, x2_bad, sigma, l),
               std::domain_error);

  EXPECT_THROW(stan::math::gp_exp_quad_cov(x_vec_1_bad, x_vec_2, sigma, l),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_exp_quad_cov(x_vec_1, x_vec_2_bad, sigma, l),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_exp_quad_cov(x_vec_1_bad, x_vec_2_bad, sigma, l),
               std::domain_error);

  EXPECT_THROW(stan::math::gp_exp_quad_cov(x_rvec_1_bad, x_rvec_2, sigma, l),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_exp_quad_cov(x_rvec_1, x_rvec_2_bad, sigma, l),
               std::domain_error);
  EXPECT_THROW(
      stan::math::gp_exp_quad_cov(x_rvec_1_bad, x_rvec_2_bad, sigma, l),
      std::domain_error);

  EXPECT_THROW(stan::math::gp_exp_quad_cov(x_vec_1_bad, x_rvec_2, sigma, l),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_exp_quad_cov(x_rvec_1_bad, x_vec_2, sigma, l),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_exp_quad_cov(x_vec_1, x_rvec_2_bad, sigma, l),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_exp_quad_cov(x_rvec_1, x_vec_2_bad, sigma, l),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_exp_quad_cov(x_vec_1_bad, x_rvec_2_bad, sigma, l),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_exp_quad_cov(x_rvec_1_bad, x_vec_2_bad, sigma, l),
               std::domain_error);
}

TEST(MathPrimMat, dim_mismatch_vec_eigen_vec_gp_exp_quad_cov2) {
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
  EXPECT_THROW(stan::math::gp_exp_quad_cov(x_vec_1, x_vec_2, sigma, l),
               std::invalid_argument);
}

TEST(MathPrimMat, vec_length_scale_eigen_gp_exp_quad_cov1) {
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
  EXPECT_NO_THROW(cov = stan::math::gp_exp_quad_cov(x, sigma, l));

  std::vector<Eigen::Matrix<double, -1, 1>> x_new(3);
  for (size_t i = 0; i < x.size(); ++i) {
    x_new[i].resize(3, 1);
    x_new[i] << 1 * i / l[0], 2 * i / l[1], 3 * i / l[2];
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_FLOAT_EQ(
          sigma * sigma
              * exp(-.5 * stan::math::squared_distance(x_new[i], x_new[j])),
          cov(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }
}

TEST(MathPrimMat, vec_length_scale_vec_eigen_vec_gp_exp_quad_cov2) {
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
  EXPECT_NO_THROW(cov = stan::math::gp_exp_quad_cov(x1, x2, sigma, l));
  EXPECT_EQ(3, cov.rows());
  EXPECT_EQ(3, cov.cols());
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_FLOAT_EQ(
          sigma * sigma
              * exp(-.5 * stan::math::squared_distance(x1_new[i], x2_new[j])),
          cov(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }
  Eigen::MatrixXd cov2;
  cov2 = stan::math::gp_exp_quad_cov(x2, x1, sigma, l);
  EXPECT_EQ(3, cov2.rows());
  EXPECT_EQ(3, cov2.cols());
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_FLOAT_EQ(
          sigma * sigma
              * exp(-.5 * stan::math::squared_distance(x2_new[i], x1_new[j])),
          cov2(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }
}

TEST(MathPrimMat, domain_error_training_sig_vec_length_scale_gp_exp_quad_cov) {
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

  double sigma_bad = -1;
  std::vector<double> l_bad(3);
  for (int i = 0; i < 3; ++i) {
    l[i] = i + 1;
  }
  l[2] = -1;

  EXPECT_THROW_MSG(stan::math::gp_exp_quad_cov(x_2, sigma, l_bad),
                   std::domain_error, " length scale");
  EXPECT_THROW_MSG(stan::math::gp_exp_quad_cov(x_2, sigma_bad, l),
                   std::domain_error, " magnitude");
  EXPECT_THROW(stan::math::gp_exp_quad_cov(x_2, sigma_bad, l_bad),
               std::domain_error);
}

TEST(MathPrimMat, nan_error_training_sig_vec_length_scale_gp_exp_quad_cov) {
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

  std::vector<Eigen::Matrix<double, -1, 1>> x_bad_2(x_2);
  x_bad_2[1](1) = std::numeric_limits<double>::quiet_NaN();

  double sigma_bad = std::numeric_limits<double>::quiet_NaN();

  std::vector<double> l_bad(3);
  for (int i = 0; i < 3; ++i) {
    l_bad[i] = i + 1;
  }
  l_bad[1] = std::numeric_limits<double>::quiet_NaN();

  EXPECT_THROW_MSG(stan::math::gp_exp_quad_cov(x_2, sigma, l_bad),
                   std::domain_error, " length scale");
  EXPECT_THROW(stan::math::gp_exp_quad_cov(x_2, sigma_bad, l_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_exp_quad_cov(x_bad_2, sigma, l_bad),
               std::domain_error);
}

TEST(MathPrimMat, nan_domain_error_gp_exp_quad_cov2_vec_length_scale) {
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

  EXPECT_THROW_MSG(stan::math::gp_exp_quad_cov(x_vec_1, x_vec_2, sigma, l_bad0),
                   std::domain_error, " length scale");
  EXPECT_THROW(stan::math::gp_exp_quad_cov(x_vec_1, x_vec_2, sigma_bad, l_bad1),
               std::domain_error);
  EXPECT_THROW_MSG(stan::math::gp_exp_quad_cov(x_vec_1, x_vec_2, sigma_bad, l),
                   std::domain_error, " magnitude");
  EXPECT_THROW_MSG(stan::math::gp_exp_quad_cov(x_vec_1, x_vec_2, sigma, l_bad1),
                   std::domain_error, " length scale");
  EXPECT_THROW_MSG(stan::math::gp_exp_quad_cov(x_vec_1, x_vec_2, sigma, l_bad2),
                   std::domain_error, " length scale");
  EXPECT_THROW_MSG(stan::math::gp_exp_quad_cov(x_vec_1, x_vec_2, sigma, l_bad3),
                   std::domain_error, " length scale");
  EXPECT_THROW_MSG(stan::math::gp_exp_quad_cov(x_vec_1, x_vec_2, sigma, l_bad4),
                   std::domain_error, " length scale");
}

TEST(MathPrimMat, zero_size_gp_exp_quad_cov) {
  double sigma = 0.2;

  std::vector<double> l(0);

  std::vector<Eigen::Matrix<double, -1, 1>> x(0);

  Eigen::MatrixXd cov;
  Eigen::MatrixXd cov2;
  cov = stan::math::gp_exp_quad_cov(x, sigma, l);
  cov2 = stan::math::gp_exp_quad_cov(x, x, sigma, l);

  EXPECT_EQ(0, cov.rows());
  EXPECT_EQ(0, cov.cols());
}

TEST(MathPrimMat, numerical_accuracy_ard_gp_exp_quad_cov) {
  double sigma = 1.0;

  std::vector<double> l(2);
  l[0] = 1;
  l[1] = 2;

  std::vector<Eigen::Matrix<double, -1, 1>> x(2);
  x[0].resize(2, 1);
  x[0] << 1, 1;
  x[1].resize(2, 1);
  x[1] << 2, 4;

  Eigen::MatrixXd cov;
  Eigen::MatrixXd cov2;
  cov = stan::math::gp_exp_quad_cov(x, sigma, l);
  cov2 = stan::math::gp_exp_quad_cov(x, x, sigma, l);

  EXPECT_FLOAT_EQ(1.0, cov(0, 0));
  EXPECT_FLOAT_EQ(1.0, cov2(0, 0));
  EXPECT_FLOAT_EQ(exp(-(1 + 9.0 / 4.0) / 2), cov(1, 0));
  EXPECT_FLOAT_EQ(exp(-(1 + 9.0 / 4.0) / 2), cov2(1, 0));
  EXPECT_FLOAT_EQ(exp(-(1 + 9.0 / 4.0) / 2), cov(0, 1));
  EXPECT_FLOAT_EQ(exp(-(1 + 9.0 / 4.0) / 2), cov2(0, 1));
  EXPECT_FLOAT_EQ(1.0, cov(1, 1));
  EXPECT_FLOAT_EQ(1.0, cov2(1, 1));
}
