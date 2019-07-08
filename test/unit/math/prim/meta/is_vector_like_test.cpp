
#include <stan/math/prim/meta.hpp>
#include <gtest/gtest.h>

#include <vector>





TEST(is_vector_like_mat, MatrixXd) {
  EXPECT_TRUE(stan::is_vector_like<Eigen::MatrixXd>::value);
}

TEST(is_vector_like_mat, vector_of_MatrixXd) {
  EXPECT_TRUE(stan::is_vector_like<std::vector<Eigen::MatrixXd> >::value);
}

TEST(is_vector_like_mat, VectorXd) {
  EXPECT_TRUE(stan::is_vector_like<Eigen::VectorXd>::value);
}

TEST(is_vector_like_mat, vector_of_VectorXd) {
  EXPECT_TRUE(stan::is_vector_like<std::vector<Eigen::VectorXd> >::value);
}

TEST(is_vector_like_mat, RowVectorXd) {
  EXPECT_TRUE(stan::is_vector_like<Eigen::RowVectorXd>::value);
}

TEST(is_vector_like_mat, vector_of_RowVectorXd) {
  EXPECT_TRUE(stan::is_vector_like<std::vector<Eigen::RowVectorXd> >::value);
}
