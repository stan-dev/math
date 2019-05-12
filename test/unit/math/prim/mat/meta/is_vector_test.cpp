#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>

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

  bool temp = is_vector<temp_matrix_d, temp_matrix_d, temp_matrix_d>::value;
  EXPECT_FALSE(temp);
  temp = is_vector<temp_matrix_d, double, double>::value;
  EXPECT_FALSE(temp);
  temp = is_vector<temp_matrix_d, temp_rowvec_d, double>::value;
  EXPECT_TRUE(temp);
  temp = is_vector<temp_matrix_d, temp_rowvec_d, temp_vec_d>::value;
  EXPECT_TRUE(temp);
  temp = is_vector<double, double, temp_vec_d>::value;
  EXPECT_TRUE(temp);
}
