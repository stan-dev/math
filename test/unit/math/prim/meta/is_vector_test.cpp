#include <stan/math/prim/meta.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MathMetaPrim, is_vector_scalar) {
  using stan::is_vector;
  EXPECT_FALSE(is_vector<double>::value);
  EXPECT_FALSE(is_vector<int>::value);
  EXPECT_FALSE(is_vector<size_t>::value);

  EXPECT_FALSE(is_vector<const double>::value);
  EXPECT_FALSE(is_vector<const int>::value);
  EXPECT_FALSE(is_vector<const size_t>::value);
  EXPECT_FALSE(is_vector<size_t*>::value);
}

TEST(MathMetaPrim, is_vector_std_vector) {
  using stan::is_vector;
  using std::vector;

  EXPECT_TRUE(is_vector<std::vector<double>>::value);
  EXPECT_TRUE(is_vector<std::vector<int>>::value);
  EXPECT_TRUE(is_vector<const std::vector<double>>::value);
  EXPECT_TRUE(is_vector<const std::vector<int>>::value);
  EXPECT_TRUE(is_vector<const std::vector<double>>::value);
  EXPECT_TRUE(is_vector<const std::vector<int>&>::value);
  EXPECT_TRUE(is_vector<std::vector<int>&&>::value);
  EXPECT_TRUE(is_vector<const std::vector<int>&&>::value);
  EXPECT_FALSE(is_vector<double>::value);
}

TEST(MathMetaPrim, is_vector_matrices) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::is_vector;

  typedef Matrix<double, Dynamic, 1> temp_vec_d;
  EXPECT_TRUE(is_vector<temp_vec_d>::value);
  EXPECT_TRUE(is_vector<const temp_vec_d>::value);

  typedef Matrix<double, 1, Dynamic> temp_rowvec_d;
  EXPECT_TRUE(is_vector<temp_rowvec_d>::value);
  EXPECT_TRUE(is_vector<const temp_rowvec_d>::value);
  EXPECT_TRUE(is_vector<const temp_rowvec_d&>::value);
  EXPECT_TRUE(is_vector<const temp_rowvec_d&&>::value);
  EXPECT_TRUE(is_vector<temp_rowvec_d&&>::value);

  typedef Matrix<double, Dynamic, Dynamic> temp_matrix_d;
  EXPECT_FALSE(is_vector<temp_matrix_d>::value);
  EXPECT_FALSE(is_vector<const temp_matrix_d>::value);
  EXPECT_FALSE(is_vector<const temp_matrix_d&>::value);
  EXPECT_FALSE(is_vector<double*>::value);
}
