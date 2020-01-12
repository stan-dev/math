#include <stan/math/prim.hpp>
#include <Eigen/Sparse>
#include <gtest/gtest.h>
#include <vector>

TEST(MathMetaPrim, primitive) {
  using stan::is_eigen;
  EXPECT_FALSE((is_eigen<bool>::value));
  EXPECT_FALSE((is_eigen<double>::value));
  EXPECT_FALSE((is_eigen<int>::value));

  EXPECT_FALSE((is_eigen<std::vector<double>>::value));

  EXPECT_TRUE((is_eigen<Eigen::EigenBase<Eigen::MatrixXd>>::value));
  EXPECT_TRUE((is_eigen<Eigen::Matrix<double, -1, -1>>::value));
  EXPECT_TRUE((is_eigen<Eigen::SparseMatrix<double>>::value));
  EXPECT_TRUE((is_eigen<Eigen::MatrixBase<Eigen::MatrixXd>>::value));

  EXPECT_TRUE((is_eigen<const Eigen::Matrix<double, -1, -1>>::value));
  EXPECT_TRUE((is_eigen<Eigen::SparseMatrix<double>&>::value));
  EXPECT_TRUE(
      (is_eigen<Eigen::MatrixBase<Eigen::Matrix<double, -1, -1>>&&>::value));

  Eigen::Matrix<double, -1, -1> a;
  Eigen::Matrix<double, -1, -1> b;

  EXPECT_TRUE((is_eigen<decltype(a * b)>::value));
  EXPECT_TRUE((is_eigen<decltype(a * b + a.transpose())>::value));
}

TEST(MathMetaPrim, expression) {
  using stan::is_eigen_matrix;
  EXPECT_FALSE((is_eigen_matrix<bool>::value));
  EXPECT_FALSE((is_eigen_matrix<double>::value));
  EXPECT_FALSE((is_eigen_matrix<int>::value));

  EXPECT_FALSE((is_eigen_matrix<std::vector<double>>::value));

  EXPECT_FALSE((is_eigen_matrix<Eigen::EigenBase<Eigen::MatrixXd>>::value));
  EXPECT_TRUE((is_eigen_matrix<Eigen::Matrix<double, -1, -1>>::value));
  Eigen::SparseMatrix<double> sparse_mat;
  EXPECT_TRUE((is_eigen_matrix<decltype(sparse_mat)>::value));
  EXPECT_FALSE((is_eigen_matrix<Eigen::MatrixBase<Eigen::MatrixXd>>::value));

  EXPECT_TRUE((is_eigen_matrix<const Eigen::Matrix<double, -1, -1>>::value));
  EXPECT_TRUE((is_eigen_matrix<Eigen::SparseMatrix<double>&>::value));
  EXPECT_FALSE((is_eigen_matrix<
                Eigen::MatrixBase<Eigen::Matrix<double, -1, -1>>&&>::value));

  Eigen::Matrix<double, -1, -1> a;
  Eigen::Matrix<double, -1, -1> b;

  EXPECT_FALSE((is_eigen_matrix<decltype(a * b)>::value));
  EXPECT_FALSE((is_eigen_matrix<decltype(a * b + a.transpose())>::value));
}
