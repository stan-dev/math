#include <stan/math/prim/meta.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(is_vector_like, double) {
  EXPECT_FALSE(stan::is_vector_like<double>::value);
}

TEST(is_vector_like, double_pointer) {
  EXPECT_TRUE(stan::is_vector_like<double *>::value);
}

TEST(is_vector_like, vector) {
  EXPECT_TRUE(stan::is_vector_like<std::vector<double>>::value);
}

TEST(is_vector_like, MatrixXd) {
  EXPECT_TRUE(stan::is_vector_like<Eigen::MatrixXd>::value);
}

TEST(is_vector_like, vector_of_MatrixXd) {
  EXPECT_TRUE(stan::is_vector_like<std::vector<Eigen::MatrixXd>>::value);
}

TEST(is_vector_like, VectorXd) {
  EXPECT_TRUE(stan::is_vector_like<Eigen::VectorXd>::value);
}

TEST(is_vector_like, vector_of_VectorXd) {
  EXPECT_TRUE(stan::is_vector_like<std::vector<Eigen::VectorXd>>::value);
}

TEST(is_vector_like, RowVectorXd) {
  EXPECT_TRUE(stan::is_vector_like<Eigen::RowVectorXd>::value);
}

TEST(is_vector_like, vector_of_RowVectorXd) {
  EXPECT_TRUE(stan::is_vector_like<std::vector<Eigen::RowVectorXd>>::value);
}

template <typename T>
class is_vector_like_test_class {
 public:
  std::vector<T> blah_{1, 2, 3, 4};
  auto operator[](T v) { return blah_[v]; }
};

TEST(is_vector_like, size_t_op) {
  EXPECT_TRUE(stan::is_vector_like<is_vector_like_test_class<size_t>>::value);
  EXPECT_TRUE(
      stan::is_vector_like<is_vector_like_test_class<unsigned int>>::value);
}
