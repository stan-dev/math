
#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <vector>







TEST(MetaTraits, is_vector) {
  using stan::is_vector;
  EXPECT_FALSE(is_vector<double>::value);
  EXPECT_FALSE(is_vector<int>::value);
  EXPECT_FALSE(is_vector<size_t>::value);

  EXPECT_FALSE(is_vector<const double>::value);
  EXPECT_FALSE(is_vector<const int>::value);
  EXPECT_FALSE(is_vector<const size_t>::value);
}




TEST(MetaTraits_arr, is_vector) {
  using stan::is_vector;
  using std::vector;

  EXPECT_TRUE(is_vector<std::vector<double> >::value);
  EXPECT_TRUE(is_vector<std::vector<int> >::value);
  EXPECT_TRUE(is_vector<std::vector<const double> >::value);
  EXPECT_TRUE(is_vector<std::vector<const int> >::value);
}



TEST(MetaTraits_mat, is_vector) {
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
