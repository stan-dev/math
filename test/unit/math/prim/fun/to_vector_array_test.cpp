#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <type_traits>
#include <vector>

TEST(ToVectorArray, values) {
  using stan::math::to_matrix;
  using stan::math::to_vector_array;

  Eigen::MatrixXd matrix(3, 2);
  matrix << 1.1, 4.4, 2.2, 5.5, 3.3, 6.6;

  auto result = to_vector_array(matrix);
  ASSERT_EQ(2, result.size());

  Eigen::VectorXd expected_col_0(3);
  expected_col_0 << 1.1, 2.2, 3.3;
  Eigen::VectorXd expected_col_1(3);
  expected_col_1 << 4.4, 5.5, 6.6;

  EXPECT_MATRIX_FLOAT_EQ(expected_col_0, result[0]);
  EXPECT_MATRIX_FLOAT_EQ(expected_col_1, result[1]);
  EXPECT_MATRIX_FLOAT_EQ(matrix, to_matrix(result));
}

TEST(ToVectorArray, emptyShapes) {
  using stan::math::to_matrix;
  using stan::math::to_vector_array;

  Eigen::MatrixXd zero_by_zero(0, 0);
  auto zero_by_zero_result = to_vector_array(zero_by_zero);
  EXPECT_EQ(0, zero_by_zero_result.size());

  Eigen::MatrixXd zero_by_three(0, 3);
  auto zero_by_three_result = to_vector_array(zero_by_three);
  ASSERT_EQ(3, zero_by_three_result.size());
  for (const auto& result : zero_by_three_result) {
    EXPECT_EQ(0, result.size());
  }
  EXPECT_MATRIX_FLOAT_EQ(zero_by_three, to_matrix(zero_by_three_result));

  Eigen::MatrixXd three_by_zero(3, 0);
  auto three_by_zero_result = to_vector_array(three_by_zero);
  EXPECT_EQ(0, three_by_zero_result.size());
}

TEST(ToVectorArray, preservesScalarType) {
  using stan::math::to_vector_array;

  Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> matrix(2, 2);
  matrix << 1, 3, 2, 4;

  auto result = to_vector_array(matrix);
  static_assert(std::is_same<decltype(result)::value_type,
                             Eigen::Matrix<int, Eigen::Dynamic, 1>>::value,
                "to_vector_array should preserve the matrix scalar type");
  ASSERT_EQ(2, result.size());
  EXPECT_EQ(1, result[0](0));
  EXPECT_EQ(2, result[0](1));
  EXPECT_EQ(3, result[1](0));
  EXPECT_EQ(4, result[1](1));
}
