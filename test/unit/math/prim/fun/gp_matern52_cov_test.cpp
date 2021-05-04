#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

TEST(MathPrimMat, vec_double_gp_matern52_cov1) {
  double sigma = 0.2;
  double l = 5.0;

  std::vector<double> x(3);
  x[0] = -2;
  x[1] = -1;
  x[2] = -0.5;

  Eigen::MatrixXd cov;
  cov = stan::math::gp_matern52_cov(x, sigma, l);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      EXPECT_FLOAT_EQ(
          sigma * sigma
              * (1
                 + std::pow(5, 0.5) / l
                       * stan::math::sqrt(
                             stan::math::squared_distance(x[i], x[j]))
                 + (5.0 / 3.0) * stan::math::squared_distance(x[i], x[j])
                       / std::pow(l, 2))
              * std::exp(
                    -1.0 * pow(5.0, 0.5)
                    * stan::math::sqrt(stan::math::squared_distance(x[i], x[j]))
                    / l),
          cov(i, j))
          << "index: (" << i << ", " << j << ")";
}

TEST(MathPrimMat, vec_eigen_gp_matern52_cov1) {
  double l = 0.2;
  double sigma = 0.3;

  std::vector<Eigen::Matrix<double, -1, 1>> x1(3);
  for (size_t i = 0; i < x1.size(); ++i) {
    x1[i].resize(3, 1);
    x1[i] << 1 * i, 2 * i, 3 * i;
  }

  Eigen::MatrixXd cov;
  EXPECT_NO_THROW(cov = stan::math::gp_matern52_cov(x1, sigma, l));
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_FLOAT_EQ(
          sigma * sigma
              * (1
                 + std::pow(5, 0.5) / l
                       * stan::math::sqrt(
                             stan::math::squared_distance(x1[i], x1[j]))
                 + (5.0 / 3.0) * stan::math::squared_distance(x1[i], x1[j])
                       / std::pow(l, 2))
              * std::exp(-1.0 * pow(5.0, 0.5)
                         * stan::math::sqrt(
                               stan::math::squared_distance(x1[i], x1[j]))
                         / l),
          cov(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }
}

TEST(MathPrimMat, vec_eigen_eigen_gp_matern52_cov1) {
  double l = 0.2;
  double sigma = 0.3;

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
  EXPECT_NO_THROW(cov = stan::math::gp_matern52_cov(x1, x2, sigma, l));
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_FLOAT_EQ(
          sigma * sigma
              * (1
                 + std::pow(5, 0.5) / l
                       * stan::math::sqrt(
                             stan::math::squared_distance(x1[i], x2[j]))
                 + (5.0 / 3.0) * stan::math::squared_distance(x1[i], x2[j])
                       / std::pow(l, 2))
              * std::exp(-1.0 * pow(5.0, 0.5)
                         * stan::math::sqrt(
                               stan::math::squared_distance(x1[i], x2[j]))
                         / l),
          cov(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }
}

TEST(MathPrimMat, vec_double_double_gp_matern52_cov1) {
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
  EXPECT_NO_THROW(cov = stan::math::gp_matern52_cov(x1, x2, sigma, l));
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      EXPECT_FLOAT_EQ(
          sigma * sigma
              * (1
                 + std::pow(5, 0.5) / l
                       * stan::math::sqrt(
                             stan::math::squared_distance(x1[i], x2[j]))
                 + (5.0 / 3.0) * stan::math::squared_distance(x1[i], x2[j])
                       / std::pow(l, 2))
              * std::exp(-1.0 * pow(5.0, 0.5)
                         * stan::math::sqrt(
                               stan::math::squared_distance(x1[i], x2[j]))
                         / l),
          cov(i, j))
          << "index: (" << i << ", " << j << ")";
}

TEST(MathPrimMat, vec_eigen_ard_gp_matern52_cov1) {
  double sigma = 0.3;

  std::vector<Eigen::Matrix<double, -1, 1>> x1(3);
  for (size_t i = 0; i < x1.size(); ++i) {
    x1[i].resize(3, 1);
    x1[i] << 1 * i, 2 * i, 3 * i;
  }

  std::vector<double> l(3);
  l[0] = 0.1;
  l[1] = 0.2;
  l[2] = 0.3;

  std::vector<Eigen::Matrix<double, -1, 1>> x_new(3);
  for (size_t i = 0; i < x_new.size(); ++i) {
    x_new[i].resize(3, 1);
    x_new[i] << 1 * i / l[0], 2 * i / l[1], 3 * i / l[2];
  }
  Eigen::MatrixXd cov;
  EXPECT_NO_THROW(cov = stan::math::gp_matern52_cov(x1, sigma, l));
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_FLOAT_EQ(
          sigma * sigma
              * (1
                 + std::pow(5.0, 0.5)
                       * stan::math::sqrt(
                             stan::math::squared_distance(x_new[i], x_new[j]))
                 + (5.0 / 3.0)
                       * stan::math::squared_distance(x_new[i], x_new[j]))
              * std::exp(-1.0 * pow(5.0, 0.5)
                         * stan::math::sqrt(stan::math::squared_distance(
                               x_new[i], x_new[j]))),
          cov(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }
}

TEST(MathPrimMat, vec_eigen_eigen_ard_gp_matern52_cov1) {
  double sigma = 0.2;

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

  std::vector<double> l(3);
  l[0] = 1;
  l[1] = 2;
  l[2] = 3;

  std::vector<Eigen::Matrix<double, -1, 1>> x1_new(3);
  for (size_t i = 0; i < x1_new.size(); ++i) {
    x1_new[i].resize(3, 1);
    x1_new[i] << 1 * i / l[0], 2 * i / l[1], 3 * i / l[2];
  }

  std::vector<Eigen::Matrix<double, -1, 1>> x2_new(3);
  for (size_t i = 0; i < x2_new.size(); ++i) {
    x2_new[i].resize(3, 1);
    x2_new[i] << 2 * i / l[0], 3 * i / l[1], 4 * i / l[2];
  }

  Eigen::MatrixXd cov;
  EXPECT_NO_THROW(cov = stan::math::gp_matern52_cov(x1, x2, sigma, l));
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_FLOAT_EQ(
          sigma * sigma
              * (1
                 + std::pow(5.0, 0.5)
                       * stan::math::sqrt(
                             stan::math::squared_distance(x1_new[i], x2_new[j]))
                 + (5.0 / 3.0)
                       * stan::math::squared_distance(x1_new[i], x2_new[j]))
              * std::exp(-1.0 * pow(5.0, 0.5)
                         * stan::math::sqrt(stan::math::squared_distance(
                               x1_new[i], x2_new[j]))),
          cov(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }
}

TEST(MathPrimMat, rvec_eigen_gp_matern52_cov1) {
  double l = 0.2;
  double sigma = 0.3;

  std::vector<Eigen::Matrix<double, -1, 1>> x1(3);
  for (size_t i = 0; i < x1.size(); ++i) {
    x1[i].resize(3, 1);
    x1[i] << 1 * i, 2 * i, 3 * i;
  }

  Eigen::MatrixXd cov;
  cov = stan::math::gp_matern52_cov(x1, sigma, l);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_FLOAT_EQ(
          sigma * sigma
              * (1
                 + std::pow(5, 0.5) / l
                       * stan::math::sqrt(
                             stan::math::squared_distance(x1[i], x1[j]))
                 + (5.0 / 3.0) * stan::math::squared_distance(x1[i], x1[j])
                       / std::pow(l, 2))
              * std::exp(-1.0 * pow(5.0, 0.5)
                         * stan::math::sqrt(
                               stan::math::squared_distance(x1[i], x1[j]))
                         / l),
          cov(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }
}

TEST(MathPrimMat, domain_err_training_sig_l_gamma_gp_matern52_cov) {
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

  EXPECT_THROW_MSG(stan::math::gp_matern32_cov(x, sigma, l_bad),
                   std::domain_error, " length scale");
  EXPECT_THROW(stan::math::gp_matern32_cov(x, sigma_bad, l_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_matern32_cov(x, sigma_bad, l_bad),
               std::domain_error);

  EXPECT_THROW_MSG(stan::math::gp_matern32_cov(x, sigma, l_bad),
                   std::domain_error, " length scale");
  EXPECT_THROW(stan::math::gp_matern32_cov(x, sigma_bad, l_bad),
               std::domain_error);

  EXPECT_THROW_MSG(stan::math::gp_matern32_cov(x_2, sigma, l_bad),
                   std::domain_error, " length scale");
  EXPECT_THROW_MSG(stan::math::gp_matern32_cov(x_2, sigma_bad, l),
                   std::domain_error, " magnitude");
  EXPECT_THROW(stan::math::gp_matern32_cov(x_2, sigma_bad, l_bad),
               std::domain_error);

  EXPECT_THROW_MSG(stan::math::gp_matern32_cov(x_2, sigma, l_bad),
                   std::domain_error, " length scale");
  EXPECT_THROW_MSG(stan::math::gp_matern32_cov(x_2, sigma_bad, l_vec),
                   std::domain_error, " magnitude");
  EXPECT_THROW(stan::math::gp_matern32_cov(x_2, sigma_bad, l_bad),
               std::domain_error);

  EXPECT_THROW_MSG(stan::math::gp_matern32_cov(x_2, sigma, l_vec_bad),
                   std::domain_error, " length scale");
  EXPECT_THROW(stan::math::gp_matern32_cov(x_2, sigma_bad, l_vec_bad),
               std::domain_error);
  EXPECT_THROW_MSG(stan::math::gp_matern32_cov(x_2, sigma_bad, l_vec),
                   std::domain_error, " magnitude");

  EXPECT_THROW_MSG(stan::math::gp_matern32_cov(x_2, x_2, sigma, l_bad),
                   std::domain_error, " length scale");
  EXPECT_THROW_MSG(stan::math::gp_matern32_cov(x_2, x_2, sigma_bad, l),
                   std::domain_error, " magnitude");
  EXPECT_THROW(stan::math::gp_matern32_cov(x_2, x_2, sigma_bad, l_bad),
               std::domain_error);

  EXPECT_THROW_MSG(stan::math::gp_matern32_cov(x_2, x_2, sigma, l_bad),
                   std::domain_error, " length scale");
  EXPECT_THROW_MSG(stan::math::gp_matern32_cov(x_2, x_2, sigma_bad, l_vec),
                   std::domain_error, " magnitude");
  EXPECT_THROW(stan::math::gp_matern32_cov(x_2, x_2, sigma_bad, l_bad),
               std::domain_error);

  EXPECT_THROW_MSG(stan::math::gp_matern32_cov(x_2, x_2, sigma, l_vec_bad),
                   std::domain_error, " length scale");
  EXPECT_THROW(stan::math::gp_matern32_cov(x_2, x_2, sigma_bad, l_vec_bad),
               std::domain_error);
  EXPECT_THROW_MSG(stan::math::gp_matern32_cov(x_2, x_2, sigma_bad, l_vec),
                   std::domain_error, " magnitude");

  EXPECT_THROW_MSG(stan::math::gp_matern32_cov(x_2, x_2, sigma, l_vec_bad),
                   std::domain_error, " length scale");
  EXPECT_THROW(stan::math::gp_matern32_cov(x_2, x_2, sigma_bad, l_vec_bad),
               std::domain_error);
  EXPECT_THROW_MSG(stan::math::gp_matern32_cov(x_2, x_2, sigma_bad, l_vec),
                   std::domain_error, " magnitude");
}

TEST(MathPrimMat, nan_error_training_sig_l_gamma_gp_matern52_cov) {
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

  std::vector<double> l_vec(2);
  l_vec[0] = 1.0;
  l_vec[1] = 2.0;
  l_vec[1] = 3.0;

  std::vector<double> l_vec_bad(2);
  l_vec_bad[0] = 1;
  l_vec_bad[1] = -1;
  l_vec_bad[1] = 2.0;

  std::vector<double> x_bad(x);
  x_bad[1] = std::numeric_limits<double>::quiet_NaN();

  std::vector<Eigen::Matrix<double, -1, 1>> x_bad_2(x_2);
  x_bad_2[1](1) = std::numeric_limits<double>::quiet_NaN();

  double sigma_bad = std::numeric_limits<double>::quiet_NaN();
  double l_bad = std::numeric_limits<double>::quiet_NaN();

  std::string msg1, msg2, msg3, msg4;
  EXPECT_THROW_MSG(stan::math::gp_matern32_cov(x, sigma, l_bad),
                   std::domain_error, " length scale");
  EXPECT_THROW_MSG(stan::math::gp_matern32_cov(x, sigma_bad, l),
                   std::domain_error, " magnitude");
  EXPECT_THROW(stan::math::gp_matern32_cov(x, sigma_bad, l_bad),
               std::domain_error);

  EXPECT_THROW(stan::math::gp_matern52_cov(x, sigma, l_bad), std::domain_error);
  EXPECT_THROW(stan::math::gp_matern52_cov(x, sigma_bad, l), std::domain_error);
  EXPECT_THROW(stan::math::gp_matern52_cov(x, sigma_bad, l_bad),
               std::domain_error);

  EXPECT_THROW(stan::math::gp_matern52_cov(x_bad, l, sigma), std::domain_error);
  EXPECT_THROW(stan::math::gp_matern52_cov(x_bad, sigma_bad, l),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_matern52_cov(x_bad, sigma_bad, l_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_matern52_cov(x_bad, sigma, l_bad),
               std::domain_error);

  EXPECT_THROW(stan::math::gp_matern52_cov(x_2, sigma, l_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_matern52_cov(x_2, sigma_bad, l),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_matern52_cov(x_2, sigma_bad, l_bad),
               std::domain_error);

  EXPECT_THROW(stan::math::gp_matern52_cov(x_bad_2, sigma, l_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_matern52_cov(x_bad_2, sigma_bad, l),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_matern52_cov(x_bad_2, sigma, l),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_matern52_cov(x_bad_2, sigma_bad, l_bad),
               std::domain_error);

  EXPECT_THROW(stan::math::gp_matern52_cov(x_2, sigma_bad, l_vec),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_matern52_cov(x_2, sigma_bad, l_vec_bad),
               std::domain_error);

  EXPECT_THROW(stan::math::gp_matern52_cov(x_bad_2, sigma, l_vec_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_matern52_cov(x_bad_2, sigma_bad, l_vec),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_matern52_cov(x_bad_2, sigma_bad, l_vec_bad),
               std::domain_error);
}

TEST(MathPrimMat, dim_mismatch_vec_eigen_rvec_gp_matern52_cov2) {
  double sigma = 0.2;
  double l = 5;

  std::vector<Eigen::Matrix<double, -1, 1>> x_vec_1(3);
  for (size_t i = 0; i < x_vec_1.size(); ++i) {
    x_vec_1[i].resize(3, 1);
    x_vec_1[i] << 1, 2, 3;
  }

  std::vector<Eigen::Matrix<double, -1, 1>> x_vec_2(3);
  for (size_t i = 0; i < x_vec_2.size(); ++i) {
    x_vec_2[i].resize(4, 1);
    x_vec_2[i] << 4, 1, 3, 1;
  }
  EXPECT_THROW(stan::math::gp_matern52_cov(x_vec_1, x_vec_2, sigma, l),
               std::invalid_argument);
}

TEST(MathPrimMat, ard_size_err_gp_matern52_cov) {
  double sigma = 0.2;

  std::vector<double> l(7);
  for (int i = 0; i < 7; ++i)
    l[i] = i + 1;

  std::vector<Eigen::Matrix<double, -1, 1>> x_vec_1(3);
  for (size_t i = 0; i < x_vec_1.size(); ++i) {
    x_vec_1[i].resize(3, 1);
    x_vec_1[i] << 1, 2, 3;
  }

  std::vector<Eigen::Matrix<double, -1, 1>> x_vec_2(3);
  for (size_t i = 0; i < x_vec_2.size(); ++i) {
    x_vec_2[i].resize(3, 1);
    x_vec_2[i] << 1, 2, 3;
  }

  EXPECT_THROW(stan::math::gp_matern52_cov(x_vec_1, sigma, l),
               std::invalid_argument);
  EXPECT_THROW(stan::math::gp_matern52_cov(x_vec_1, x_vec_2, sigma, l),
               std::invalid_argument);
}

TEST(MathPrimMat, zero_size_gp_matern52_cov) {
  double sigma = 0.2;

  std::vector<double> l(0);

  std::vector<Eigen::Matrix<double, -1, 1>> x(0);

  Eigen::MatrixXd cov;
  Eigen::MatrixXd cov2;
  EXPECT_NO_THROW(cov = stan::math::gp_matern52_cov(x, sigma, l));
  EXPECT_NO_THROW(cov2 = stan::math::gp_matern52_cov(x, x, sigma, l));
  EXPECT_EQ(0, cov.rows());
  EXPECT_EQ(0, cov.cols());
}
