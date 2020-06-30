#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <vector>
#include <stdexcept>

TEST(MathMatrix, to_array_1d_matrix) {
  using stan::math::to_array_1d;
  Eigen::MatrixXd a(3, 3);
  a << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  std::vector<double> a_correct{1, 4, 7, 2, 5, 8, 3, 6, 9};
  std::vector<double> a_res = to_array_1d(a);
  EXPECT_STD_VECTOR_FLOAT_EQ(a_res, a_correct);
}

TEST(MathMatrix, to_array_1d_matrix_block) {
  using stan::math::to_array_1d;
  Eigen::MatrixXd a(3, 3);
  a << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  std::vector<double> a_correct{2, 5, 3, 6};
  std::vector<double> a_res = to_array_1d(a.block(0, 1, 2, 2));
  EXPECT_STD_VECTOR_FLOAT_EQ(a_res, a_correct);
}
