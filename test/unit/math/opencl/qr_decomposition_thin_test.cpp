#ifdef STAN_OPENCL

#include <stan/math/opencl/prim.hpp>
#include <gtest/gtest.h>
#include <algorithm>

using stan::math::matrix_cl;

TEST(MathMatrixGPU, qr_decomposition_thin_prim_small) {
  double eps = 1e-10;
  int rows = 4;
  int cols = 3;
  Eigen::MatrixXd m(rows, cols);
  m << 2, 2, 2, 2, 3, 4, 5, 6, -1, -2, 0, 3;

  matrix_cl<double> m_cl(m);

  matrix_cl<double> q_cl = stan::math::qr_thin_Q(m_cl);
  matrix_cl<double> r_cl = stan::math::qr_thin_R(m_cl);

  Eigen::MatrixXd q = stan::math::from_matrix_cl(q_cl);
  Eigen::MatrixXd r = stan::math::from_matrix_cl(r_cl);

  EXPECT_EQ(r.rows(), std::min(rows, cols));
  EXPECT_EQ(r.cols(), cols);
  EXPECT_EQ(q.rows(), rows);
  EXPECT_EQ(q.cols(), std::min(rows, cols));

  Eigen::MatrixXd reconstructed_m = q * r;
  Eigen::MatrixXd expected_identity = q.transpose() * q;

  Eigen::MatrixXd identity = Eigen::MatrixXd::Identity(rows, cols);

  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      EXPECT_NEAR(reconstructed_m(i, j), m(i, j), eps);
      if (i > j && i < std::min(rows, cols)) {
        EXPECT_NEAR(r(i, j), 0, eps);
      }
    }
    if (i < std::min(rows, cols)) {
      for (int j = 0; j < std::min(rows, cols); j++) {
        EXPECT_NEAR(expected_identity(i, j), identity(i, j), eps);
      }
    }
  }
}

TEST(MathMatrixGPU, qr_decomposition_thin_prim_big) {
  double eps = 1e-10;
  int rows_opts[] = {1, 1, 317, 317, 411};
  int cols_opts[] = {1, 317, 1, 411, 317};
  for (int k = 0; k < 5; k++) {
    int rows = rows_opts[k];
    int cols = cols_opts[k];
    Eigen::MatrixXd m = Eigen::MatrixXd::Random(rows, cols);

    matrix_cl<double> m_cl(m);

    matrix_cl<double> q_cl = stan::math::qr_thin_Q(m_cl);
    matrix_cl<double> r_cl = stan::math::qr_thin_R(m_cl);

    Eigen::MatrixXd q = stan::math::from_matrix_cl(q_cl);
    Eigen::MatrixXd r = stan::math::from_matrix_cl(r_cl);

    EXPECT_EQ(r.rows(), std::min(rows, cols));
    EXPECT_EQ(r.cols(), cols);
    EXPECT_EQ(q.rows(), rows);
    EXPECT_EQ(q.cols(), std::min(rows, cols));

    Eigen::MatrixXd reconstructed_m = q * r;
    Eigen::MatrixXd expected_identity = q.transpose() * q;

    Eigen::MatrixXd identity = Eigen::MatrixXd::Identity(rows, rows);

    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < cols; j++) {
        EXPECT_NEAR(reconstructed_m(i, j), m(i, j), eps);
        if (i > j && i < std::min(rows, cols)) {
          EXPECT_NEAR(r(i, j), 0, eps);
        }
      }
      if (i < std::min(rows, cols)) {
        for (int j = 0; j < std::min(rows, cols); j++) {
          EXPECT_NEAR(expected_identity(i, j), identity(i, j), eps);
        }
      }
    }
  }
}

#endif
