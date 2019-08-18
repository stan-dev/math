#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MetaTraits, is_vector) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::is_vector;

  typedef Matrix<double, Dynamic, 1> temp_vec_d;
  EXPECT_TRUE(is_vector<temp_vec_d>::value);
  EXPECT_TRUE(is_vector<const temp_vec_d>::value);

  typedef Matrix<double, 1, Dynamic> temp_rowvec_d;
  EXPECT_TRUE(is_vector<temp_rowvec_d>::value);
  EXPECT_TRUE(is_vector<const temp_rowvec_d>::value);

  typedef Matrix<double, Dynamic, Dynamic> temp_matrix_d;
  EXPECT_FALSE(is_vector<temp_matrix_d>::value);
  EXPECT_FALSE(is_vector<const temp_matrix_d>::value);
}


TEST(MetaTraits, is_eigen_vector_vector) {
  EXPECT_FALSE(stan::is_eigen_row_vector<std::vector<double>>::value);
  EXPECT_FALSE(stan::is_eigen_col_vector<std::vector<double>>::value);
  EXPECT_FALSE(stan::is_eigen_vector<std::vector<double>>::value);

  EXPECT_FALSE(stan::is_eigen_row_vector<const std::vector<double>>::value);
  EXPECT_FALSE(stan::is_eigen_col_vector<const std::vector<double>>::value);
  EXPECT_FALSE(stan::is_eigen_vector<const std::vector<double>>::value);
}

TEST(MetaTraits, is_eigen_vector_eigen_row_vector) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  typedef Eigen::Matrix<double, 1, Dynamic> rowvec_d;
  EXPECT_TRUE(stan::is_eigen_row_vector<rowvec_d>::value);
  EXPECT_TRUE(stan::is_eigen_row_vector<const rowvec_d>::value);
  EXPECT_FALSE(stan::is_eigen_col_vector<rowvec_d>::value);
  EXPECT_FALSE(stan::is_eigen_col_vector<const rowvec_d>::value);
  EXPECT_TRUE(stan::is_eigen_vector<rowvec_d>::value);
  EXPECT_TRUE(stan::is_eigen_vector<const rowvec_d>::value);
}

TEST(MetaTraits, is_eigen_vector_eigen_col_vector) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  typedef Eigen::Matrix<double, Dynamic, 1> vec_d;
  EXPECT_FALSE(stan::is_eigen_row_vector<vec_d>::value);
  EXPECT_FALSE(stan::is_eigen_row_vector<const vec_d>::value);
  EXPECT_TRUE(stan::is_eigen_col_vector<vec_d>::value);
  EXPECT_TRUE(stan::is_eigen_col_vector<const vec_d>::value);
  EXPECT_TRUE(stan::is_eigen_vector<vec_d>::value);
  EXPECT_TRUE(stan::is_eigen_vector<const vec_d>::value);
}

TEST(MetaTraits, is_eigen_vector_eigen_matrix_vector) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  typedef Eigen::Matrix<double, Dynamic, Dynamic> matrix_d;
  EXPECT_FALSE(stan::is_eigen_row_vector<matrix_d>::value);
  EXPECT_FALSE(stan::is_eigen_row_vector<const matrix_d>::value);
  EXPECT_FALSE(stan::is_eigen_col_vector<matrix_d>::value);
  EXPECT_FALSE(stan::is_eigen_col_vector<const matrix_d>::value);
  EXPECT_FALSE(stan::is_eigen_vector<matrix_d>::value);
  EXPECT_FALSE(stan::is_eigen_vector<const matrix_d>::value);
}
