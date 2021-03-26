#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(MathMatrixPrimMat, rep_matrix) {
  using stan::math::rep_matrix;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> x
      = rep_matrix(2.0, 3, 4);
  EXPECT_EQ(12, x.size());
  EXPECT_EQ(3, x.rows());
  EXPECT_EQ(4, x.cols());
  for (int i = 0; i < x.rows(); ++i)
    for (int j = 0; j < x.cols(); ++j)
      EXPECT_FLOAT_EQ(2.0, x(i, j));

  EXPECT_THROW(rep_matrix(2.0, -1, 3), std::domain_error);
  EXPECT_THROW(rep_matrix(2.0, 3, -1), std::domain_error);
}

TEST(MathMatrixPrimMat, rep_matrix_always_returns_double) {
  int a = 1;
  auto x = stan::math::rep_matrix(a, 3, 4);

  EXPECT_TRUE((std::is_same<double, stan::value_type_t<decltype(x)>>::value));
}

TEST(MathMatrixPrimMat, rep_matrix_vec) {
  using stan::math::rep_matrix;
  Eigen::Matrix<double, Eigen::Dynamic, 1> v(3);
  v << 1.0, 4.0, 9.0;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> x = rep_matrix(v, 4);
  EXPECT_EQ(12, x.size());
  EXPECT_EQ(3, x.rows());
  EXPECT_EQ(4, x.cols());
  for (int i = 0; i < x.rows(); ++i)
    for (int j = 0; j < x.cols(); ++j)
      EXPECT_FLOAT_EQ(v(i), x(i, j));

  EXPECT_THROW(rep_matrix(v, -1), std::domain_error);
}

TEST(MathMatrixPrimMat, rep_matrix_row_vec) {
  using stan::math::rep_matrix;
  Eigen::Matrix<double, 1, Eigen::Dynamic> rv(3);
  rv << 1.0, 4.0, 9.0;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> x = rep_matrix(rv, 4);
  EXPECT_EQ(12, x.size());
  EXPECT_EQ(4, x.rows());
  EXPECT_EQ(3, x.cols());
  for (int i = 0; i < x.rows(); ++i)
    for (int j = 0; j < x.cols(); ++j)
      EXPECT_FLOAT_EQ(rv(j), x(i, j));

  EXPECT_THROW(rep_matrix(rv, -1), std::domain_error);
}
