#include <gtest/gtest.h>
#include <stan/math/prim/mat.hpp>
#include <limits>
#include <string>
#include <vector>

TEST(MathPrimMat, 1d_one_x) {
  double sigma_squared = 0.25;

  std::vector<double> x(3);
  x[0] = -2;
  x[1] = -1;
  x[2] = -0.5;

  Eigen::MatrixXd ref(3, 3);

  ref << 1.00, 0.500, 0.2500, 0.50, 0.250, 0.1250, 0.25, 0.125, 0.0625;

  Eigen::MatrixXd cov;
  EXPECT_NO_THROW(cov = stan::math::gp_dot_prod_cov(x, sigma_squared));

  for (int i = 0; i < ref.rows(); i++)
    for (int j = 0; j < ref.cols(); j++)
      EXPECT_FLOAT_EQ(ref(i, j), cov(i, j))
          << "index: (" << i << ", " << j << ")";
}

TEST(MathPrimMat, 1d_two_x) {
  double sigma_squared = 0.25;

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
  ref << -0.500, -1.00, -1.500, -2.0, -0.250, -0.50, -0.750, -1.0, -0.125,
      -0.25, -0.375, -0.5;

  Eigen::MatrixXd cov;
  EXPECT_NO_THROW(cov = stan::math::gp_dot_prod_cov(x1, x2, sigma_squared));

  EXPECT_EQ(ref.rows(), cov.rows());
  EXPECT_EQ(ref.cols(), cov.cols());
  for (int i = 0; i < ref.rows(); i++)
    for (int j = 0; j < ref.cols(); j++)
      EXPECT_FLOAT_EQ(ref(i, j), cov(i, j))
          << "index: (" << i << ", " << j << ")";
}

TEST(MathPrimMat, 2d_scalar_one_x) {
  double sigma_squared = 1.1;

  Eigen::MatrixXd ref(4, 4);

  ref << 14.3, 28.6, 42.9, 57.2,
    28.6, 57.2, 85.8, 114.4,
    42.9, 85.8, 128.7, 171.6,
    57.2, 114.4, 171.6, 228.8;

  std::vector<Eigen::Matrix<double, -1, 1>> x(4);
  for (size_t i = 0; i < x.size(); ++i) {
    x[i].resize(2, 1);
    x[i] << 2 * (i + 1), 3 * (i + 1);
  }

  Eigen::MatrixXd cov;
  EXPECT_NO_THROW(cov = stan::math::gp_dot_prod_cov(x, sigma_squared));
  for (int i = 0; i < ref.rows(); i++) {
    for (int j = 0; j < ref.cols(); j++) {
      EXPECT_FLOAT_EQ(ref(i, j), cov(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }
}

TEST(MathPrimMat, 2d_scalar_two_x) {
  double sigma_squared = 1.1;

  Eigen::MatrixXd ref(3, 4);

  ref << -15.4, -23.1, -30.8, -38.5,
    -30.8, -46.2, -61.6, -77.0,
    -46.2, -69.3, -92.4, -115.5;

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
  EXPECT_NO_THROW(cov = stan::math::gp_dot_prod_cov(x1, x2, sigma_squared));
  EXPECT_EQ(ref.rows(), cov.rows());
  EXPECT_EQ(ref.cols(), cov.cols());
  for (int i = 0; i < ref.rows(); ++i) {
    for (int j = 0; j < ref.cols(); ++j) {
      EXPECT_FLOAT_EQ(ref(i, j), cov(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }
}

TEST(MathPrimMat, 2d_vector_one_x) {
  Eigen::VectorXd diagonal_Sigma(2);

  diagonal_Sigma << 1.5, 1.7;

  Eigen::MatrixXd ref(4, 4);

  ref << 21.3, 42.6, 63.9, 85.2,
    42.6, 85.2, 127.8, 170.4,
    63.9, 127.8, 191.7, 255.6,
    85.2, 170.4, 255.6, 340.8;

  std::vector<Eigen::Matrix<double, -1, 1>> x(4);
  for (size_t i = 0; i < x.size(); ++i) {
    x[i].resize(2, 1);
    x[i] << 2 * (i + 1), 3 * (i + 1);
  }

  Eigen::MatrixXd cov;
  EXPECT_NO_THROW(cov = stan::math::gp_dot_prod_cov(x, diagonal_Sigma));
  for (int i = 0; i < ref.rows(); i++) {
    for (int j = 0; j < ref.cols(); j++) {
      EXPECT_FLOAT_EQ(ref(i, j), cov(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }
}

TEST(MathPrimMat, 2d_vector_two_x) {
  Eigen::VectorXd diagonal_Sigma(3);

  diagonal_Sigma << 1.1, 2.70, 4.70;

  Eigen::MatrixXd ref(3, 4);

  ref << -68.2, -133.5, -198.8, -264.1,
    -136.4, -267.0, -397.6, -528.2,
    -204.6, -400.5, -596.4, -792.3;

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
  EXPECT_NO_THROW(cov = stan::math::gp_dot_prod_cov(x1, x2, diagonal_Sigma));
  EXPECT_EQ(ref.rows(), cov.rows());
  EXPECT_EQ(ref.cols(), cov.cols());
  for (int i = 0; i < ref.rows(); ++i) {
    for (int j = 0; j < ref.cols(); ++j) {
      EXPECT_FLOAT_EQ(ref(i, j), cov(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }
}

TEST(MathPrimMat, 2d_matrix_one_x) {
  Eigen::MatrixXd Sigma(2, 2);

  Sigma << 1.5, 0.5, 0.5, 1.7;

  Eigen::MatrixXd ref(4, 4);

  ref << 27.3, 54.6, 81.9, 109.2, 54.6, 109.2, 163.8, 218.4, 81.9, 163.8, 245.7,
      327.6, 109.2, 218.4, 327.6, 436.8;

  std::vector<Eigen::Matrix<double, -1, 1>> x(4);
  for (size_t i = 0; i < x.size(); ++i) {
    x[i].resize(2, 1);
    x[i] << 2 * (i + 1), 3 * (i + 1);
  }

  Eigen::MatrixXd cov;
  EXPECT_NO_THROW(cov = stan::math::gp_dot_prod_cov(x, Sigma));
  for (int i = 0; i < ref.rows(); i++) {
    for (int j = 0; j < ref.cols(); j++) {
      EXPECT_FLOAT_EQ(ref(i, j), cov(i, j))
          << "index: (" << i << ", " << j << ")";
    }
  }
}

TEST(MathPrimMat, 2d_matrix_two_x) {
  Eigen::MatrixXd Sigma(3, 3);

  Sigma << 1.1, 0.30, 0.10, 0.3, 2.70, 0.25, 0.1, 0.25, 4.70;

  Eigen::MatrixXd ref(3, 4);

  ref << -70.45, -133.1, -195.75, -258.4, -140.90, -266.2, -391.50, -516.8,
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

TEST(MathPrimMat, 1d_one_x_error) {
  double sigma_squared = 0.5;
  double sigma_squared_nan = std::numeric_limits<double>::quiet_NaN();
  double sigma_squared_inf = std::numeric_limits<double>::infinity();
  double sigma_squared_not_positive = 0.0;
  double sigma_squared_negative = -1.0;

  std::vector<double> x = { -2, -1, -0.5 };
  std::vector<double> x_nan = { -2, -1, std::numeric_limits<double>::quiet_NaN() };
  std::vector<double> x_inf = { -2, -1, std::numeric_limits<double>::infinity() };

  Eigen::MatrixXd cov;
  EXPECT_THROW(cov = stan::math::gp_dot_prod_cov(x_nan, sigma_squared), std::domain_error);
  EXPECT_THROW(cov = stan::math::gp_dot_prod_cov(x_inf, sigma_squared), std::domain_error);
  EXPECT_THROW(cov = stan::math::gp_dot_prod_cov(x, sigma_squared_nan), std::domain_error);
  EXPECT_THROW(cov = stan::math::gp_dot_prod_cov(x, sigma_squared_inf), std::domain_error);
  EXPECT_THROW(cov = stan::math::gp_dot_prod_cov(x, sigma_squared_not_positive), std::domain_error);
  EXPECT_THROW(cov = stan::math::gp_dot_prod_cov(x, sigma_squared_negative), std::domain_error);
}

TEST(MathPrimMat, 1d_two_x_error) {
  double sigma_squared = 0.5;
  double sigma_squared_nan = std::numeric_limits<double>::quiet_NaN();
  double sigma_squared_inf = std::numeric_limits<double>::infinity();
  double sigma_squared_not_positive = 0.0;
  double sigma_squared_negative = -1.0;

  std::vector<double> x = { -2, -1, -0.5 };
  std::vector<double> x_nan = { -2, -1, std::numeric_limits<double>::quiet_NaN() };
  std::vector<double> x_inf = { -2, -1, std::numeric_limits<double>::infinity() };

  Eigen::MatrixXd cov;
  EXPECT_THROW(cov = stan::math::gp_dot_prod_cov(x_nan, x, sigma_squared), std::domain_error);
  EXPECT_THROW(cov = stan::math::gp_dot_prod_cov(x_inf, x, sigma_squared), std::domain_error);
  EXPECT_THROW(cov = stan::math::gp_dot_prod_cov(x, x_nan, sigma_squared), std::domain_error);
  EXPECT_THROW(cov = stan::math::gp_dot_prod_cov(x, x_inf, sigma_squared), std::domain_error);
  EXPECT_THROW(cov = stan::math::gp_dot_prod_cov(x, x, sigma_squared_nan), std::domain_error);
  EXPECT_THROW(cov = stan::math::gp_dot_prod_cov(x, x, sigma_squared_inf), std::domain_error);
  EXPECT_THROW(cov = stan::math::gp_dot_prod_cov(x, x, sigma_squared_not_positive), std::domain_error);
  EXPECT_THROW(cov = stan::math::gp_dot_prod_cov(x, x, sigma_squared_negative), std::domain_error);
}

TEST(MathPrimMat, 2d_scalar_one_x_error) {
  double sigma_squared = 0.5;
  double sigma_squared_nan = std::numeric_limits<double>::quiet_NaN();
  double sigma_squared_inf = std::numeric_limits<double>::infinity();
  double sigma_squared_not_positive = 0.0;
  double sigma_squared_negative = -1.0;

  std::vector<Eigen::Matrix<double, -1, 1>> x(1);
  x[0].resize(2, 1);
  x[0] << 1, 2;
  std::vector<Eigen::Matrix<double, -1, 1>> x_nan(1);
  x_nan[0].resize(2, 1);
  x_nan[0] << 1, std::numeric_limits<double>::quiet_NaN();
  std::vector<Eigen::Matrix<double, -1, 1>> x_inf(1);
  x_inf[0].resize(2, 1);
  x_inf[0] << 1, std::numeric_limits<double>::infinity();

  Eigen::MatrixXd cov;
  EXPECT_THROW(cov = stan::math::gp_dot_prod_cov(x_nan, sigma_squared), std::domain_error);
  EXPECT_THROW(cov = stan::math::gp_dot_prod_cov(x_inf, sigma_squared), std::domain_error);
  EXPECT_THROW(cov = stan::math::gp_dot_prod_cov(x, sigma_squared_nan), std::domain_error);
  EXPECT_THROW(cov = stan::math::gp_dot_prod_cov(x, sigma_squared_inf), std::domain_error);
  EXPECT_THROW(cov = stan::math::gp_dot_prod_cov(x, sigma_squared_not_positive), std::domain_error);
  EXPECT_THROW(cov = stan::math::gp_dot_prod_cov(x, sigma_squared_negative), std::domain_error);
}

TEST(MathPrimMat, 2d_scalar_two_x_error) {
  double sigma_squared = 0.5;
  double sigma_squared_nan = std::numeric_limits<double>::quiet_NaN();
  double sigma_squared_inf = std::numeric_limits<double>::infinity();
  double sigma_squared_not_positive = 0.0;
  double sigma_squared_negative = -1.0;

  std::vector<Eigen::Matrix<double, -1, 1>> x(1);
  x[0].resize(2, 1);
  x[0] << 1, 2;
  std::vector<Eigen::Matrix<double, -1, 1>> x_nan(1);
  x_nan[0].resize(2, 1);
  x_nan[0] << 1, std::numeric_limits<double>::quiet_NaN();
  std::vector<Eigen::Matrix<double, -1, 1>> x_inf(1);
  x_inf[0].resize(2, 1);
  x_inf[0] << 1, std::numeric_limits<double>::infinity();
  std::vector<Eigen::Matrix<double, -1, 1>> x_size(1);
  x_size[0].resize(3, 1);
  x_size[0] << 1, 2, 3;

  Eigen::MatrixXd cov;
  EXPECT_THROW(cov = stan::math::gp_dot_prod_cov(x_nan, x, sigma_squared), std::domain_error);
  EXPECT_THROW(cov = stan::math::gp_dot_prod_cov(x_inf, x, sigma_squared), std::domain_error);
  EXPECT_THROW(cov = stan::math::gp_dot_prod_cov(x_size, x, sigma_squared), std::invalid_argument);
  EXPECT_THROW(cov = stan::math::gp_dot_prod_cov(x, x_nan, sigma_squared), std::domain_error);
  EXPECT_THROW(cov = stan::math::gp_dot_prod_cov(x, x_inf, sigma_squared), std::domain_error);
  EXPECT_THROW(cov = stan::math::gp_dot_prod_cov(x, x_size, sigma_squared), std::invalid_argument);
  EXPECT_THROW(cov = stan::math::gp_dot_prod_cov(x, x, sigma_squared_nan), std::domain_error);
  EXPECT_THROW(cov = stan::math::gp_dot_prod_cov(x, x, sigma_squared_inf), std::domain_error);
  EXPECT_THROW(cov = stan::math::gp_dot_prod_cov(x, x, sigma_squared_not_positive), std::domain_error);
  EXPECT_THROW(cov = stan::math::gp_dot_prod_cov(x, x, sigma_squared_negative), std::domain_error);
}

TEST(MathPrimMat, 2d_vector_one_x_error) {
  Eigen::VectorXd diagonal_Sigma(2);
  diagonal_Sigma << 1.5, 1.7;
  Eigen::VectorXd diagonal_Sigma_nan(2);
  diagonal_Sigma_nan << 1.5, std::numeric_limits<double>::quiet_NaN();
  Eigen::VectorXd diagonal_Sigma_inf(2);
  diagonal_Sigma_inf << 1.5, std::numeric_limits<double>::infinity();
  Eigen::VectorXd diagonal_Sigma_not_positive(2);
  diagonal_Sigma_not_positive << 1.5, 0.0;
  Eigen::VectorXd diagonal_Sigma_negative(2);
  diagonal_Sigma_negative << 1.5, -1.0;
  Eigen::VectorXd diagonal_Sigma_size(3);
  diagonal_Sigma_size << 1.5, 1.7, 1.1;

  std::vector<Eigen::Matrix<double, -1, 1>> x(1);
  x[0].resize(2, 1);
  x[0] << 1, 2;
  std::vector<Eigen::Matrix<double, -1, 1>> x_nan(1);
  x_nan[0].resize(2, 1);
  x_nan[0] << 1, std::numeric_limits<double>::quiet_NaN();
  std::vector<Eigen::Matrix<double, -1, 1>> x_inf(1);
  x_inf[0].resize(2, 1);
  x_inf[0] << 1, std::numeric_limits<double>::infinity();
  std::vector<Eigen::Matrix<double, -1, 1>> x_size(1);
  x_size[0].resize(3, 1);
  x_size[0] << 1, 2, 3;

  Eigen::MatrixXd cov;
  EXPECT_THROW(cov = stan::math::gp_dot_prod_cov(x_nan, diagonal_Sigma), std::domain_error);
  EXPECT_THROW(cov = stan::math::gp_dot_prod_cov(x_inf, diagonal_Sigma), std::domain_error);
  EXPECT_THROW(cov = stan::math::gp_dot_prod_cov(x_size, diagonal_Sigma), std::invalid_argument);
  EXPECT_THROW(cov = stan::math::gp_dot_prod_cov(x, diagonal_Sigma_nan), std::domain_error);
  EXPECT_THROW(cov = stan::math::gp_dot_prod_cov(x, diagonal_Sigma_inf), std::domain_error);
  EXPECT_THROW(cov = stan::math::gp_dot_prod_cov(x, diagonal_Sigma_not_positive), std::domain_error);
  EXPECT_THROW(cov = stan::math::gp_dot_prod_cov(x, diagonal_Sigma_negative), std::domain_error);
  EXPECT_THROW(cov = stan::math::gp_dot_prod_cov(x, diagonal_Sigma_size), std::invalid_argument);
}

TEST(MathPrimMat, 2d_vector_two_x_error) {
  Eigen::VectorXd diagonal_Sigma(2);
  diagonal_Sigma << 1.5, 1.7;
  Eigen::VectorXd diagonal_Sigma_nan(2);
  diagonal_Sigma_nan << 1.5, std::numeric_limits<double>::quiet_NaN();
  Eigen::VectorXd diagonal_Sigma_inf(2);
  diagonal_Sigma_inf << 1.5, std::numeric_limits<double>::infinity();
  Eigen::VectorXd diagonal_Sigma_not_positive(2);
  diagonal_Sigma_not_positive << 1.5, 0.0;
  Eigen::VectorXd diagonal_Sigma_negative(2);
  diagonal_Sigma_negative << 1.5, -1.0;
  Eigen::VectorXd diagonal_Sigma_size(3);
  diagonal_Sigma_size << 1.5, 1.7, 1.1;

  std::vector<Eigen::Matrix<double, -1, 1>> x(1);
  x[0].resize(2, 1);
  x[0] << 1, 2;
  std::vector<Eigen::Matrix<double, -1, 1>> x_nan(1);
  x_nan[0].resize(2, 1);
  x_nan[0] << 1, std::numeric_limits<double>::quiet_NaN();
  std::vector<Eigen::Matrix<double, -1, 1>> x_inf(1);
  x_inf[0].resize(2, 1);
  x_inf[0] << 1, std::numeric_limits<double>::infinity();
  std::vector<Eigen::Matrix<double, -1, 1>> x_size(1);
  x_size[0].resize(3, 1);
  x_size[0] << 1, 2, 3;

  Eigen::MatrixXd cov;
  EXPECT_THROW(cov = stan::math::gp_dot_prod_cov(x_nan, x, diagonal_Sigma), std::domain_error);
  EXPECT_THROW(cov = stan::math::gp_dot_prod_cov(x_inf, x, diagonal_Sigma), std::domain_error);
  EXPECT_THROW(cov = stan::math::gp_dot_prod_cov(x_size, x, diagonal_Sigma), std::invalid_argument);
  EXPECT_THROW(cov = stan::math::gp_dot_prod_cov(x, x_nan, diagonal_Sigma), std::domain_error);
  EXPECT_THROW(cov = stan::math::gp_dot_prod_cov(x, x_inf, diagonal_Sigma), std::domain_error);
  EXPECT_THROW(cov = stan::math::gp_dot_prod_cov(x, x_size, diagonal_Sigma), std::invalid_argument);
  EXPECT_THROW(cov = stan::math::gp_dot_prod_cov(x, x, diagonal_Sigma_nan), std::domain_error);
  EXPECT_THROW(cov = stan::math::gp_dot_prod_cov(x, x, diagonal_Sigma_inf), std::domain_error);
  EXPECT_THROW(cov = stan::math::gp_dot_prod_cov(x, x, diagonal_Sigma_not_positive), std::domain_error);
  EXPECT_THROW(cov = stan::math::gp_dot_prod_cov(x, x, diagonal_Sigma_negative), std::domain_error);
  EXPECT_THROW(cov = stan::math::gp_dot_prod_cov(x, x, diagonal_Sigma_size), std::invalid_argument);
}

TEST(MathPrimMat, 2d_matrix_one_x_error) {
  Eigen::MatrixXd Sigma(2, 2);
  Sigma << 1.5, 0.1, 0.1, 1.7;
  Eigen::MatrixXd Sigma_nan(2, 2);
  Sigma_nan << 1.5, 0.1, 0.1, std::numeric_limits<double>::quiet_NaN();
  Eigen::MatrixXd Sigma_inf(2, 2);
  Sigma_inf << 1.5, 0.1, 0.1, std::numeric_limits<double>::infinity();
  Eigen::MatrixXd Sigma_not_pos_def(2, 2);
  Sigma_not_pos_def << 1.5, 0.0, 0.0, -1.0;
  Eigen::MatrixXd Sigma_size(3, 3);
  Sigma_size << 1.5, 0.0, 0.0, 0.0, 1.7, 0.0, 0.0, 0.0, 1.1;

  std::vector<Eigen::Matrix<double, -1, 1>> x(1);
  x[0].resize(2, 1);
  x[0] << 1, 2;
  std::vector<Eigen::Matrix<double, -1, 1>> x_nan(1);
  x_nan[0].resize(2, 1);
  x_nan[0] << 1, std::numeric_limits<double>::quiet_NaN();
  std::vector<Eigen::Matrix<double, -1, 1>> x_inf(1);
  x_inf[0].resize(2, 1);
  x_inf[0] << 1, std::numeric_limits<double>::infinity();
  std::vector<Eigen::Matrix<double, -1, 1>> x_size(1);
  x_size[0].resize(3, 1);
  x_size[0] << 1, 2, 3;

  Eigen::MatrixXd cov;
  EXPECT_THROW(cov = stan::math::gp_dot_prod_cov(x_nan, Sigma), std::domain_error);
  EXPECT_THROW(cov = stan::math::gp_dot_prod_cov(x_inf, Sigma), std::domain_error);
  EXPECT_THROW(cov = stan::math::gp_dot_prod_cov(x_size, Sigma), std::invalid_argument);
  EXPECT_THROW(cov = stan::math::gp_dot_prod_cov(x, Sigma_nan), std::domain_error);
  EXPECT_THROW(cov = stan::math::gp_dot_prod_cov(x, Sigma_inf), std::domain_error);
  EXPECT_THROW(cov = stan::math::gp_dot_prod_cov(x, Sigma_not_pos_def), std::domain_error);
  EXPECT_THROW(cov = stan::math::gp_dot_prod_cov(x, Sigma_size), std::invalid_argument);
}

TEST(MathPrimMat, 2d_matrix_two_x_error) {
  Eigen::MatrixXd Sigma(2, 2);
  Sigma << 1.5, 0.1, 0.1, 1.7;
  Eigen::MatrixXd Sigma_nan(2, 2);
  Sigma_nan << 1.5, 0.1, 0.1, std::numeric_limits<double>::quiet_NaN();
  Eigen::MatrixXd Sigma_inf(2, 2);
  Sigma_inf << 1.5, 0.1, 0.1, std::numeric_limits<double>::infinity();
  Eigen::MatrixXd Sigma_not_pos_def(2, 2);
  Sigma_not_pos_def << 1.5, 0.0, 0.0, -1.0;
  Eigen::MatrixXd Sigma_size(3, 3);
  Sigma_size << 1.5, 0.0, 0.0, 0.0, 1.7, 0.0, 0.0, 0.0, 1.1;

  std::vector<Eigen::Matrix<double, -1, 1>> x(1);
  x[0].resize(2, 1);
  x[0] << 1, 2;
  std::vector<Eigen::Matrix<double, -1, 1>> x_nan(1);
  x_nan[0].resize(2, 1);
  x_nan[0] << 1, std::numeric_limits<double>::quiet_NaN();
  std::vector<Eigen::Matrix<double, -1, 1>> x_inf(1);
  x_inf[0].resize(2, 1);
  x_inf[0] << 1, std::numeric_limits<double>::infinity();
  std::vector<Eigen::Matrix<double, -1, 1>> x_size(1);
  x_size[0].resize(3, 1);
  x_size[0] << 1, 2, 3;

  Eigen::MatrixXd cov;
  EXPECT_THROW(cov = stan::math::gp_dot_prod_cov(x_nan, x, Sigma), std::domain_error);
  EXPECT_THROW(cov = stan::math::gp_dot_prod_cov(x_inf, x, Sigma), std::domain_error);
  EXPECT_THROW(cov = stan::math::gp_dot_prod_cov(x_size, x, Sigma), std::invalid_argument);
  EXPECT_THROW(cov = stan::math::gp_dot_prod_cov(x, x_nan, Sigma), std::domain_error);
  EXPECT_THROW(cov = stan::math::gp_dot_prod_cov(x, x_inf, Sigma), std::domain_error);
  EXPECT_THROW(cov = stan::math::gp_dot_prod_cov(x, x_size, Sigma), std::invalid_argument);
  EXPECT_THROW(cov = stan::math::gp_dot_prod_cov(x, x, Sigma_nan), std::domain_error);
  EXPECT_THROW(cov = stan::math::gp_dot_prod_cov(x, x, Sigma_inf), std::domain_error);
  EXPECT_THROW(cov = stan::math::gp_dot_prod_cov(x, x, Sigma_not_pos_def), std::domain_error);
  EXPECT_THROW(cov = stan::math::gp_dot_prod_cov(x, x, Sigma_size), std::invalid_argument);
}
