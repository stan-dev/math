#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(is_eigen, MatrixXd) {
  EXPECT_TRUE(stan::is_eigen<Eigen::MatrixXd>::value);
}

TEST(is_eigen, vector_of_MatrixXd) {
  EXPECT_FALSE(stan::is_eigen<std::vector<Eigen::MatrixXd>>::value);
}

TEST(is_eigen, VectorXd) {
  EXPECT_TRUE(stan::is_eigen<Eigen::VectorXd>::value);
}

TEST(is_eigen, vector_of_VectorXd) {
  EXPECT_FALSE(stan::is_eigen<std::vector<Eigen::VectorXd>>::value);
}

TEST(is_eigen, RowVectorXd) {
  EXPECT_TRUE(stan::is_eigen<Eigen::RowVectorXd>::value);
}

TEST(is_eigen, vector_of_RowVectorXd) {
  EXPECT_FALSE(stan::is_eigen<std::vector<Eigen::RowVectorXd>>::value);
}

TEST(is_eigen, arithmetic_type) { EXPECT_FALSE(stan::is_eigen<double>::value); }
