#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>

TEST(MetaTraits, contains_vector) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::contains_vector;

  typedef Matrix<double, Dynamic, 1> temp_vec_d;
  EXPECT_TRUE(contains_vector<temp_vec_d>::value);
  EXPECT_TRUE(contains_vector<const temp_vec_d>::value);

  typedef Matrix<double, 1, Dynamic> temp_rowvec_d;
  EXPECT_TRUE(contains_vector<temp_rowvec_d>::value);
  EXPECT_TRUE(contains_vector<const temp_rowvec_d>::value);

  typedef Matrix<double, Dynamic, Dynamic> temp_matrix_d;
  EXPECT_FALSE(contains_vector<temp_matrix_d>::value);
  EXPECT_FALSE(contains_vector<const temp_matrix_d>::value);

  bool temp
      = contains_vector<temp_vec_d, temp_vec_d, double, temp_matrix_d>::value;
  EXPECT_TRUE(temp);

  temp = contains_vector<double, temp_matrix_d, temp_matrix_d>::value;
  EXPECT_FALSE(temp);

  temp = contains_vector<double, double, temp_matrix_d, double, temp_matrix_d,
                         double, double>::value;
  EXPECT_FALSE(temp);

  temp = contains_vector<double, double, temp_matrix_d, double, temp_matrix_d,
                         double, temp_vec_d>::value;
  EXPECT_TRUE(temp);
}
