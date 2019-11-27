#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/scal/meta/is_vector.hpp>
#include <Eigen/Sparse>
#include <gtest/gtest.h>
#include <vector>

TEST(MathMeta, is_std_vector) {
  using stan::is_std_vector;
  EXPECT_FALSE((is_std_vector<bool>::value));
  EXPECT_FALSE((is_std_vector<double>::value));
  EXPECT_FALSE((is_std_vector<int>::value));

  EXPECT_FALSE((is_std_vector<std::vector<double>>::value));
  EXPECT_FALSE((is_std_vector<const std::vector<double>>::value));
  EXPECT_FALSE((is_std_vector<const std::vector<double>&>::value));

  EXPECT_FALSE((is_std_vector<Eigen::EigenBase<Eigen::MatrixXd>>::value));
  EXPECT_FALSE((is_std_vector<Eigen::Matrix<double, -1, -1>>::value));
  EXPECT_FALSE((is_std_vector<Eigen::SparseMatrix<double>>::value));
  EXPECT_FALSE((is_std_vector<Eigen::MatrixBase<Eigen::MatrixXd>>::value));

  EXPECT_FALSE((is_std_vector<const Eigen::Matrix<double, -1, -1>>::value));
  EXPECT_FALSE((is_std_vector<Eigen::SparseMatrix<double>&>::value));
  EXPECT_FALSE((is_std_vector<
                Eigen::MatrixBase<Eigen::Matrix<double, -1, -1>>&&>::value));

  Eigen::Matrix<double, -1, -1> a;
  Eigen::Matrix<double, -1, -1> b;

  EXPECT_FALSE((is_std_vector<decltype(a * b)>::value));
  EXPECT_FALSE((is_std_vector<decltype(a * b + a.transpose())>::value));
}
