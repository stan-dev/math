
#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <vector>







TEST(MetaTraits, contains_vector) {
  using stan::contains_vector;
  EXPECT_FALSE(contains_vector<double>::value);
  EXPECT_FALSE(contains_vector<int>::value);
  EXPECT_FALSE(contains_vector<size_t>::value);

  EXPECT_FALSE(contains_vector<const double>::value);
  EXPECT_FALSE(contains_vector<const int>::value);
  EXPECT_FALSE(contains_vector<const size_t>::value);
}




TEST(MetaTraits_arr, contains_vector) {
  using stan::contains_vector;
  using std::vector;

  EXPECT_TRUE(contains_vector<std::vector<double> >::value);
  EXPECT_TRUE(contains_vector<std::vector<int> >::value);
  EXPECT_TRUE(contains_vector<std::vector<const double> >::value);
  EXPECT_TRUE(contains_vector<std::vector<const int> >::value);
}



TEST(MetaTraits_mat, contains_vector) {
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
