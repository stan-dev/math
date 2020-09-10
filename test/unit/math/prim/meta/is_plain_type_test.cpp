#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MathMetaPrim, is_plain_type) {
  using stan::is_plain_type;
  EXPECT_TRUE((is_plain_type<bool>::value));
  EXPECT_TRUE((is_plain_type<double>::value));
  EXPECT_TRUE((is_plain_type<int>::value));

  EXPECT_TRUE((is_plain_type<std::vector<double>>::value));
  EXPECT_TRUE((is_plain_type<const std::vector<double>>::value));
  EXPECT_TRUE((is_plain_type<const std::vector<double>&>::value));

  EXPECT_TRUE((is_plain_type<Eigen::Matrix<double, -1, -1>>::value));
  EXPECT_TRUE((is_plain_type<Eigen::SparseMatrix<double>>::value));

  EXPECT_TRUE((is_plain_type<const Eigen::Matrix<double, -1, -1>>::value));
  EXPECT_TRUE((is_plain_type<Eigen::SparseMatrix<double>&>::value));
  Eigen::Matrix<double, -1, -1> a;
  Eigen::Matrix<double, -1, -1> b;

  EXPECT_FALSE((is_plain_type<decltype(a * b)>::value));
  EXPECT_FALSE((is_plain_type<decltype(a * b + a.transpose())>::value));
  EXPECT_FALSE((is_plain_type<
                Eigen::MatrixBase<Eigen::Matrix<double, -1, -1>>&&>::value));
  EXPECT_FALSE((is_plain_type<Eigen::MatrixBase<Eigen::MatrixXd>>::value));
}
