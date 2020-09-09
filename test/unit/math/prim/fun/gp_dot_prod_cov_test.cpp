#include <gtest/gtest.h>
#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <limits>
#include <string>
#include <vector>

TEST(MathPrimMat, vec_double_gp_dot_prod_cov0) {
  double sigma = 0.5;
  double sigma_sq = pow(sigma, 2);

  std::vector<double> x(3);
  x[0] = -2;
  x[1] = -1;
  x[2] = -0.5;

  Eigen::MatrixXd cov;
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

TEST(MathPrimMat, vec_NaN_x_gp_dot_prod_cov_cov0) {
  double sigma = 1;
  std::vector<double> x(2);
  x[0] = 1;
  x[1] = std::numeric_limits<double>::quiet_NaN();

  EXPECT_THROW_MSG(stan::math::gp_dot_prod_cov(x, sigma), std::domain_error,
                   " x");
}

TEST(MathPrimMat, vec_NaN_sigma_gp_dot_prod_cov_cov0) {
  double sigma = std::numeric_limits<double>::quiet_NaN();
  std::vector<double> x(2);
  x[0] = 1;
  x[1] = 2;

  EXPECT_THROW_MSG(stan::math::gp_dot_prod_cov(x, sigma), std::domain_error,
                   " sigma");
}

TEST(MathPrimMat, vec_NaN_x1_gp_dot_prod_cov_cov0) {
  double sigma = 2;

  std::vector<double> x1(2);
  x1[0] = std::numeric_limits<double>::quiet_NaN();
  x1[1] = 2;

  std::vector<double> x2(2);
  x2[0] = 1;
  x2[1] = 2;

  EXPECT_THROW_MSG(stan::math::gp_dot_prod_cov(x1, x2, sigma),
                   std::domain_error, " x1");
}

TEST(MathPrimMat, vec_NaN_x2_gp_dot_prod_cov_cov0) {
  double sigma = 2;

  std::vector<double> x1(2);
  x1[0] = 1;
  x1[1] = 2;

  std::vector<double> x2(2);
  x2[0] = 1;
  x2[1] = std::numeric_limits<double>::quiet_NaN();

  EXPECT_THROW_MSG(stan::math::gp_dot_prod_cov(x1, x2, sigma),
                   std::domain_error, " x2");
}

TEST(MathPrimMat, vec_NaN_x1_x2_sigma_gp_dot_prod_cov_cov0) {
  double sigma = std::numeric_limits<double>::quiet_NaN();

  std::vector<double> x1(2);
  x1[0] = 1;
  x1[1] = 2;

  std::vector<double> x2(2);
  x2[0] = 1;
  x2[1] = 3;

  EXPECT_THROW_MSG(stan::math::gp_dot_prod_cov(x1, x2, sigma),
                   std::domain_error, " sigma");
}

TEST(MathPrimMat, vec_inf_x_gp_dot_prod_cov_cov0) {
  double sigma = 1;
  std::vector<double> x(2);
  x[0] = 1;
  x[1] = std::numeric_limits<double>::infinity();

  EXPECT_THROW_MSG(stan::math::gp_dot_prod_cov(x, sigma), std::domain_error,
                   " x");
}

TEST(MathPrimMat, vec_inf_sigma_gp_dot_prod_cov_cov0) {
  double sigma = std::numeric_limits<double>::infinity();
  std::vector<double> x(2);
  x[0] = 1;
  x[1] = 2;

  EXPECT_THROW_MSG(stan::math::gp_dot_prod_cov(x, sigma), std::domain_error,
                   " sigma");
}

TEST(MathPrimMat, vec_inf_x1_gp_dot_prod_cov_cov0) {
  double sigma = 2;

  std::vector<double> x1(2);
  x1[0] = std::numeric_limits<double>::infinity();
  x1[1] = 2;

  std::vector<double> x2(2);
  x2[0] = 1;
  x2[1] = 2;

  EXPECT_THROW_MSG(stan::math::gp_dot_prod_cov(x1, x2, sigma),
                   std::domain_error, " x1");
}

TEST(MathPrimMat, vec_inf_x2_gp_dot_prod_cov_cov0) {
  double sigma = 2;

  std::vector<double> x1(2);
  x1[0] = 1;
  x1[1] = 2;

  std::vector<double> x2(2);
  x2[0] = 1;
  x2[1] = std::numeric_limits<double>::infinity();

  EXPECT_THROW_MSG(stan::math::gp_dot_prod_cov(x1, x2, sigma),
                   std::domain_error, " x2");
}

TEST(MathPrimMat, vec_inf_x1_x2_sigma_gp_dot_prod_cov_cov0) {
  double sigma = std::numeric_limits<double>::infinity();

  std::vector<double> x1(2);
  x1[0] = 1;
  x1[1] = 2;

  std::vector<double> x2(2);
  x2[0] = 1;
  x2[1] = 3;

  EXPECT_THROW_MSG(stan::math::gp_dot_prod_cov(x1, x2, sigma),
                   std::domain_error, " sigma");
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
