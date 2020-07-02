#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

TEST(MathPrimMat, vec_double_gp_exponential_cov1) {
  double sigma = 0.2;
  double l = 5.0;

  std::vector<double> x = {-2, -1, -0.5};

  Eigen::MatrixXd cov;
  cov = stan::math::gp_exponential_cov(x, sigma, l);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      EXPECT_FLOAT_EQ(
          sigma * sigma * std::exp(-stan::math::distance(x[i], x[j]) / (l)),
          cov(i, j))
          << "index: (" << i << ", " << j << ")";
}

TEST(MathPrimMat, vec_double_double_gp_exponential_cov1) {
  double sigma = 0.2;
  double l = 5.0;

  std::vector<double> x1 = {-2, -1, -0.5};

  std::vector<double> x2 = {3, -4, -0.7};

  Eigen::MatrixXd cov;
  cov = stan::math::gp_exponential_cov(x1, x2, sigma, l);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      EXPECT_FLOAT_EQ(
          sigma * sigma
              * std::exp(-1.0 * stan::math::distance(x1[i], x2[j]) / l),
          cov(i, j))
          << "index: (" << i << ", " << j << ")";
}

TEST(MathPrimMat, vec_eigen_gp_exponential_cov1) {
  double l = 0.2;
  double sigma = 0.3;

  std::vector<Eigen::Matrix<double, -1, 1>> x1(3);
  for (size_t i = 0; i < x1.size(); ++i) {
    x1[i].resize(3, 1);
    x1[i] << 1 * i, 2 * i, 3 * i;
  }

  Eigen::MatrixXd cov;
  cov = stan::math::gp_exponential_cov(x1, sigma, l);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_FLOAT_EQ(
          sigma * sigma
              * std::exp(-1.0 * stan::math::distance(x1[i], x1[j]) / l),
          cov(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }
}

TEST(MathPrimMat, vec_eigen_vec_eigen_gp_exponential_cov1) {
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
  cov = stan::math::gp_exponential_cov(x1, x2, sigma, l);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_FLOAT_EQ(
          sigma * sigma
              * std::exp(-1.0
                         * stan::math::sqrt(
                               stan::math::squared_distance(x1[i], x2[j]))
                         / l),
          cov(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }

  cov = stan::math::gp_exponential_cov(x2, x1, sigma, l);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_FLOAT_EQ(
          sigma * sigma
              * std::exp(-1.0 * stan::math::distance(x2[i], x1[j]) / l),
          cov(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }
}

TEST(MathPrimMat, vec_eigen_vec_eigen_ard_gp_exponential_cov1) {
  double sigma = 0.3;
  double temp;

  std::vector<double> l = {1, 2, 3};

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
  cov = stan::math::gp_exponential_cov(x1, x2, sigma, l);
  std::vector<Eigen::Matrix<double, -1, 1>> x1_new
      = stan::math::divide_columns(x1, l);
  std::vector<Eigen::Matrix<double, -1, 1>> x2_new
      = stan::math::divide_columns(x2, l);

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_FLOAT_EQ(
          sigma * sigma
              * std::exp(-1.0 * stan::math::distance(x1_new[i], x2_new[j])),
          cov(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }
  cov = stan::math::gp_exponential_cov(x2, x1, sigma, l);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_FLOAT_EQ(
          sigma * sigma * std::exp(-stan::math::distance(x2_new[i], x1_new[j])),
          cov(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }
}

TEST(MathPrimMat, domain_err_training_sig_l_gp_exp_cov) {
  double sigma = .2;
  double l = 7;
  double sigma_bad = -1.0;
  double l_bad = -1.0;

  std::vector<double> l_vec = {1, 2, 3};

  std::vector<double> l_vec_bad = {1, -1, 2};

  std::vector<double> x = {-2, -1, -0.5};

  std::vector<Eigen::Matrix<double, -1, 1>> x_2(3);
  for (size_t i = 0; i < x_2.size(); ++i) {
    x_2[i].resize(3, 1);
    x_2[i] << 1, 2, 3;
  }

  EXPECT_THROW_MSG(stan::math::gp_exponential_cov(x, sigma, l_bad),
                   std::domain_error, " length scale");

  EXPECT_THROW(stan::math::gp_exponential_cov(x, sigma_bad, l_bad),
               std::domain_error);

  EXPECT_THROW_MSG(stan::math::gp_exponential_cov(x, sigma, l_bad),
                   std::domain_error, " length scale");
  EXPECT_THROW(stan::math::gp_exponential_cov(x, sigma_bad, l_bad),
               std::domain_error);

  EXPECT_THROW_MSG(stan::math::gp_exponential_cov(x_2, sigma, l_bad),
                   std::domain_error, " length scale");
  EXPECT_THROW(stan::math::gp_exponential_cov(x_2, sigma_bad, l_bad),
               std::domain_error);

  EXPECT_THROW_MSG(stan::math::gp_exponential_cov(x_2, sigma, l_vec_bad),
                   std::domain_error, " length scale");
  EXPECT_THROW(stan::math::gp_exponential_cov(x_2, sigma_bad, l_vec_bad),
               std::domain_error);

  EXPECT_THROW_MSG(stan::math::gp_exponential_cov(x_2, sigma_bad, l_vec),
                   std::domain_error, " magnitude");

  EXPECT_THROW_MSG(stan::math::gp_exponential_cov(x_2, x_2, sigma, l_bad),
                   std::domain_error, " length scale");
  EXPECT_THROW_MSG(stan::math::gp_exponential_cov(x_2, x_2, sigma_bad, l),
                   std::domain_error, " magnitude");
  EXPECT_THROW(stan::math::gp_exponential_cov(x_2, x_2, sigma_bad, l_bad),
               std::domain_error);

  EXPECT_THROW_MSG(stan::math::gp_exponential_cov(x_2, x_2, sigma, l_bad),
                   std::domain_error, " length scale");
  EXPECT_THROW_MSG(stan::math::gp_exponential_cov(x_2, x_2, sigma_bad, l_vec),
                   std::domain_error, " magnitude");
  EXPECT_THROW(stan::math::gp_exponential_cov(x_2, x_2, sigma_bad, l_bad),
               std::domain_error);

  EXPECT_THROW_MSG(stan::math::gp_exponential_cov(x_2, x_2, sigma, l_vec_bad),
                   std::domain_error, " length scale");
  EXPECT_THROW(stan::math::gp_exponential_cov(x_2, x_2, sigma_bad, l_vec_bad),
               std::domain_error);
  EXPECT_THROW_MSG(stan::math::gp_exponential_cov(x_2, x_2, sigma_bad, l_vec),
                   std::domain_error, " magnitude");

  EXPECT_THROW_MSG(stan::math::gp_exponential_cov(x_2, x_2, sigma, l_vec_bad),
                   std::domain_error, " length scale");
  EXPECT_THROW(stan::math::gp_exponential_cov(x_2, x_2, sigma_bad, l_vec_bad),
               std::domain_error);
  EXPECT_THROW_MSG(stan::math::gp_exponential_cov(x_2, x_2, sigma_bad, l_vec),
                   std::domain_error, " magnitude");

  EXPECT_THROW_MSG(stan::math::gp_exponential_cov(x_2, x_2, sigma, l_vec_bad),
                   std::domain_error, " length scale");
  EXPECT_THROW(stan::math::gp_exponential_cov(x_2, x_2, sigma_bad, l_vec_bad),
               std::domain_error);
  EXPECT_THROW_MSG(stan::math::gp_exponential_cov(x_2, x_2, sigma_bad, l_vec),
                   std::domain_error, " magnitude");
}

TEST(MathPrimMat, nan_error_training_sig_l_gp_exp_cov) {
  double sigma = 0.2;
  double l = 5;

  std::vector<double> x = {-2, -1, -0.5};

  std::vector<Eigen::Matrix<double, -1, 1>> x_2(3);
  for (size_t i = 0; i < x_2.size(); ++i) {
    x_2[i].resize(3, 1);
    x_2[i] << 1, 2, 3;
  }

  std::vector<double> l_vec = {1, 2, 3};

  std::vector<double> l_vec_bad = {1, -1, 3};

  std::vector<double> x_bad(x);
  x_bad[1] = std::numeric_limits<double>::quiet_NaN();

  std::vector<Eigen::Matrix<double, -1, 1>> x_bad_2(x_2);
  x_bad_2[1](1) = std::numeric_limits<double>::quiet_NaN();

  double sigma_bad = std::numeric_limits<double>::quiet_NaN();
  double l_bad = std::numeric_limits<double>::quiet_NaN();

  EXPECT_THROW_MSG(stan::math::gp_exponential_cov(x, sigma, l_bad),
                   std::domain_error, " length scale");
  EXPECT_THROW_MSG(stan::math::gp_exponential_cov(x, sigma_bad, l),
                   std::domain_error, " magnitude");
  EXPECT_THROW(stan::math::gp_exponential_cov(x, sigma_bad, l_bad),
               std::domain_error);

  EXPECT_THROW(stan::math::gp_exponential_cov(x_bad, l, sigma),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_exponential_cov(x_bad, sigma_bad, l),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_exponential_cov(x_bad, sigma_bad, l_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_exponential_cov(x_bad, sigma, l_bad),
               std::domain_error);

  EXPECT_THROW(stan::math::gp_exponential_cov(x_2, sigma, l_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_exponential_cov(x_2, sigma_bad, l),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_exponential_cov(x_2, sigma_bad, l_bad),
               std::domain_error);

  EXPECT_THROW(stan::math::gp_exponential_cov(x_bad_2, sigma, l_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_exponential_cov(x_bad_2, sigma_bad, l),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_exponential_cov(x_bad_2, sigma, l),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_exponential_cov(x_bad_2, sigma_bad, l_bad),
               std::domain_error);

  EXPECT_THROW(stan::math::gp_exponential_cov(x_2, sigma, l_vec_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_exponential_cov(x_2, sigma_bad, l_vec),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_exponential_cov(x_2, sigma_bad, l_vec_bad),
               std::domain_error);

  EXPECT_THROW(stan::math::gp_exponential_cov(x_bad_2, sigma, l_vec_bad),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_exponential_cov(x_bad_2, sigma_bad, l_vec),
               std::domain_error);
  EXPECT_THROW(stan::math::gp_exponential_cov(x_bad_2, sigma_bad, l_vec_bad),
               std::domain_error);
}

TEST(MathPrimMat, check_dim_mismatch_gp_exp_cov) {
  double sig = 1.0;
  double l = 1.0;

  std::vector<Eigen::Matrix<double, -1, 1>> x(2);
  x[0].resize(2, 1);
  x[0] << 1, 2;
  x[1].resize(3, 1);
  x[1] << 1, 2, 3;

  EXPECT_THROW(stan::math::gp_exponential_cov(x, sig, l),
               std::invalid_argument);

  std::vector<Eigen::Matrix<double, -1, 1>> x1(2);
  x1[0].resize(2, 1);
  x1[0] << 1, 2;
  x1[1].resize(2, 1);
  x1[1] << 1, 2;

  std::vector<Eigen::Matrix<double, -1, 1>> x2(3);
  x2[0].resize(2, 1);
  x2[0] << 1, 2;
  x2[1].resize(3, 1);
  x2[1] << 1, 2, 3;

  EXPECT_THROW(stan::math::gp_exponential_cov(x1, x2, sig, l),
               std::invalid_argument);
}

TEST(MathPrimMat, zero_size_gp_exp_cov) {
  double sigma = 0.2;

  std::vector<double> l(0);

  std::vector<Eigen::Matrix<double, -1, 1>> x(0);

  Eigen::MatrixXd cov;
  Eigen::MatrixXd cov2;
  EXPECT_NO_THROW(cov = stan::math::gp_exponential_cov(x, sigma, l));
  EXPECT_NO_THROW(cov2 = stan::math::gp_exponential_cov(x, x, sigma, l));
  EXPECT_EQ(0, cov.rows());
  EXPECT_EQ(0, cov.cols());
}

TEST(MathPrimMat, calculations_gp_exp_cov) {
  double sigma = 1.0;
  double l = 1.0;

  std::vector<double> x1 = {1, 2};
  std::vector<double> x2 = {2, 3};

  Eigen::MatrixXd cov;
  EXPECT_NO_THROW(cov = stan::math::gp_exponential_cov(x1, sigma, l));
  ASSERT_FLOAT_EQ(1.0, cov(0, 0));
  ASSERT_FLOAT_EQ(1.0, cov(1, 1));
  ASSERT_FLOAT_EQ(exp(-1), cov(1, 0));
  ASSERT_FLOAT_EQ(exp(-1), cov(0, 1));
  EXPECT_NO_THROW(cov = stan::math::gp_exponential_cov(x1, x2, sigma, l));
  ASSERT_FLOAT_EQ(exp(-1), cov(0, 0));
  ASSERT_FLOAT_EQ(exp(-1), cov(1, 1));
  ASSERT_FLOAT_EQ(1.0, cov(1, 0));
  ASSERT_FLOAT_EQ(exp(-2), cov(0, 1));
}

TEST(MathPrimMat, calculations_ard_gp_exp_cov) {
  double sigma = 1.0;

  std::vector<double> l = {1, 2};

  std::vector<Eigen::Matrix<double, -1, 1>> x(2);
  x[0].resize(2, 1);
  x[0] << 1, 1;
  x[1].resize(2, 1);
  x[1] << 2, 4;

  Eigen::MatrixXd cov;
  Eigen::MatrixXd cov2;
  cov = stan::math::gp_exponential_cov(x, sigma, l);
  cov2 = stan::math::gp_exponential_cov(x, x, sigma, l);

  EXPECT_FLOAT_EQ(1.0, cov(0, 0));
  EXPECT_FLOAT_EQ(1.0, cov2(0, 0));
  EXPECT_FLOAT_EQ(exp(-sqrt(1 + 9.0 / 4.0)), cov(1, 0));
  EXPECT_FLOAT_EQ(exp(-sqrt(1 + 9.0 / 4.0)), cov2(1, 0));
  EXPECT_FLOAT_EQ(exp(-sqrt(1 + 9.0 / 4.0)), cov(0, 1));
  EXPECT_FLOAT_EQ(exp(-sqrt(1 + 9.0 / 4.0)), cov2(0, 1));
  EXPECT_FLOAT_EQ(1.0, cov(1, 1));
  EXPECT_FLOAT_EQ(1.0, cov2(1, 1));
}
