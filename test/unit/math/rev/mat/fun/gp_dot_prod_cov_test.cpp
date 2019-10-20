#include <gtest/gtest.h>
#include <stan/math/rev/mat.hpp>
#include <test/unit/math/rev/mat/util.hpp>
#include <limits>
#include <string>
#include <vector>

TEST(RevMath, 1d_x) {
  double sigma_squared = 0.25;
  stan::math::var sigma_squaredv = 0.25;

  std::vector<double> x = { -2, -1, -0.5 };
  std::vector<stan::math::var> xv = { -2, -1, -0.5 };

  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> cov;
  EXPECT_NO_THROW(cov = stan::math::gp_dot_prod_cov(x, sigma_squaredv));
  EXPECT_NO_THROW(cov = stan::math::gp_dot_prod_cov(xv, sigma_squared));
  EXPECT_NO_THROW(cov = stan::math::gp_dot_prod_cov(xv, sigma_squaredv));

  EXPECT_NO_THROW(cov = stan::math::gp_dot_prod_cov(x, x, sigma_squaredv));
  EXPECT_NO_THROW(cov = stan::math::gp_dot_prod_cov(x, xv, sigma_squared));
  EXPECT_NO_THROW(cov = stan::math::gp_dot_prod_cov(x, xv, sigma_squaredv));
  EXPECT_NO_THROW(cov = stan::math::gp_dot_prod_cov(xv, x, sigma_squaredv));
  EXPECT_NO_THROW(cov = stan::math::gp_dot_prod_cov(xv, xv, sigma_squared));
  EXPECT_NO_THROW(cov = stan::math::gp_dot_prod_cov(xv, xv, sigma_squaredv));
}

TEST(RevMath, nd_scalar_x) {
  double sigma_squared = 0.25;
  stan::math::var sigma_squaredv = 0.25;

  std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1>> x(4);
  std::vector<Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1>> xv(x.size());
  for (size_t i = 0; i < x.size(); ++i) {
    x[i].resize(2, 1);
    xv[i].resize(2, 1);

    x[i] << 2 * (i + 1), 3 * (i + 1);
    xv[i] << 2 * (i + 1), 3 * (i + 1);
  }


  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> cov;
  EXPECT_NO_THROW(cov = stan::math::gp_dot_prod_cov(x, sigma_squaredv));
  EXPECT_NO_THROW(cov = stan::math::gp_dot_prod_cov(xv, sigma_squared));
  EXPECT_NO_THROW(cov = stan::math::gp_dot_prod_cov(xv, sigma_squaredv));

  EXPECT_NO_THROW(cov = stan::math::gp_dot_prod_cov(x, x, sigma_squaredv));
  EXPECT_NO_THROW(cov = stan::math::gp_dot_prod_cov(x, xv, sigma_squared));
  EXPECT_NO_THROW(cov = stan::math::gp_dot_prod_cov(x, xv, sigma_squaredv));
  EXPECT_NO_THROW(cov = stan::math::gp_dot_prod_cov(xv, x, sigma_squaredv));
  EXPECT_NO_THROW(cov = stan::math::gp_dot_prod_cov(xv, xv, sigma_squared));
  EXPECT_NO_THROW(cov = stan::math::gp_dot_prod_cov(xv, xv, sigma_squaredv));
}

TEST(RevMath, nd_vector_x) {
  Eigen::VectorXd diagonal_Sigma(2);
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> diagonal_Sigmav(2);
  diagonal_Sigma << 1.5, 1.7;
  diagonal_Sigmav << 1.5, 1.7;

  std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1>> x(4);
  std::vector<Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1>> xv(x.size());
  for (size_t i = 0; i < x.size(); ++i) {
    x[i].resize(2, 1);
    xv[i].resize(2, 1);

    x[i] << 2 * (i + 1), 3 * (i + 1);
    xv[i] << 2 * (i + 1), 3 * (i + 1);
  }


  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> cov;
  EXPECT_NO_THROW(cov = stan::math::gp_dot_prod_cov(x, diagonal_Sigmav));
  EXPECT_NO_THROW(cov = stan::math::gp_dot_prod_cov(xv, diagonal_Sigma));
  EXPECT_NO_THROW(cov = stan::math::gp_dot_prod_cov(xv, diagonal_Sigmav));

  EXPECT_NO_THROW(cov = stan::math::gp_dot_prod_cov(x, x, diagonal_Sigmav));
  EXPECT_NO_THROW(cov = stan::math::gp_dot_prod_cov(x, xv, diagonal_Sigma));
  EXPECT_NO_THROW(cov = stan::math::gp_dot_prod_cov(x, xv, diagonal_Sigmav));
  EXPECT_NO_THROW(cov = stan::math::gp_dot_prod_cov(xv, x, diagonal_Sigmav));
  EXPECT_NO_THROW(cov = stan::math::gp_dot_prod_cov(xv, xv, diagonal_Sigma));
  EXPECT_NO_THROW(cov = stan::math::gp_dot_prod_cov(xv, xv, diagonal_Sigmav));
}


TEST(RevMath, nd_matrix_x) {
  Eigen::MatrixXd Sigma(2, 2);
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> Sigmav(2, 2);
  Sigma << 1.5, 0.1, 0.1, 1.7;
  Sigmav << 1.5, 0.1, 0.1, 1.7;

  std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1>> x(4);
  std::vector<Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1>> xv(x.size());
  for (size_t i = 0; i < x.size(); ++i) {
    x[i].resize(2, 1);
    xv[i].resize(2, 1);

    x[i] << 2 * (i + 1), 3 * (i + 1);
    xv[i] << 2 * (i + 1), 3 * (i + 1);
  }


  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> cov;
  EXPECT_NO_THROW(cov = stan::math::gp_dot_prod_cov(x, Sigmav));
  EXPECT_NO_THROW(cov = stan::math::gp_dot_prod_cov(xv, Sigma));
  EXPECT_NO_THROW(cov = stan::math::gp_dot_prod_cov(xv, Sigmav));

  EXPECT_NO_THROW(cov = stan::math::gp_dot_prod_cov(x, x, Sigmav));
  EXPECT_NO_THROW(cov = stan::math::gp_dot_prod_cov(x, xv, Sigma));
  EXPECT_NO_THROW(cov = stan::math::gp_dot_prod_cov(x, xv, Sigmav));
  EXPECT_NO_THROW(cov = stan::math::gp_dot_prod_cov(xv, x, Sigmav));
  EXPECT_NO_THROW(cov = stan::math::gp_dot_prod_cov(xv, xv, Sigma));
  EXPECT_NO_THROW(cov = stan::math::gp_dot_prod_cov(xv, xv, Sigmav));
}

