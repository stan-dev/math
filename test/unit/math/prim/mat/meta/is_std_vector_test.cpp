#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <Eigen/Sparse>
#include <gtest/gtest.h>
#include <vector>

TEST(MathMeta, is_std_vector) {
  using stan::is_std_vector;
  EXPECT_FALSE((is_std_vector<bool>::value));
  EXPECT_FALSE((is_std_vector<double>::value));
  EXPECT_FALSE((is_std_vector<int>::value));

  EXPECT_TRUE((is_std_vector<std::vector<double>>::value));
  EXPECT_TRUE((is_std_vector<const std::vector<double>>::value));
  EXPECT_TRUE((is_std_vector<const std::vector<double>&>::value));

  EXPECT_FALSE((is_std_vector<Eigen::EigenBase<Eigen::MatrixXd>>::value));
  EXPECT_FALSE((is_std_vector<Eigen::Matrix<double, -1, -1>>::value));
  EXPECT_FALSE((is_std_vector<Eigen::SparseMatrix<double>>::value));
  EXPECT_FALSE((is_std_vector<Eigen::MatrixBase<Eigen::MatrixXd>>::value));

  EXPECT_FALSE((is_std_vector<const Eigen::Matrix<double, -1, -1>>::value));
  EXPECT_FALSE((is_std_vector<Eigen::SparseMatrix<double>&>::value));
  EXPECT_FALSE((is_std_vector<
                Eigen::MatrixBase<Eigen::Matrix<double, -1, -1>>&&>::value));
}
