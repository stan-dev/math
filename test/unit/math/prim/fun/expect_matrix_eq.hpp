#ifndef TEST_UNIT_MATH_PRIM_FUN_EXPECT_MATRIX_EQ_HPP
#define TEST_UNIT_MATH_PRIM_FUN_EXPECT_MATRIX_EQ_HPP

#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <vector>

void expect_matrix_eq(
    const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& a,
    const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& b) {
  EXPECT_EQ(a.rows(), b.rows());
  EXPECT_EQ(a.cols(), b.cols());
  for (int i = 0; i < a.rows(); ++i)
    for (int j = 0; j < a.cols(); ++j)
      EXPECT_FLOAT_EQ(a(i, j), b(i, j));
}

template <typename T, typename = stan::require_arithmetic_t<T>>
void expect_std_vector_eq(const std::vector<T>& a, const std::vector<T>& b) {
  EXPECT_EQ(a.size(), b.size());
  for (int i = 0; i < a.size(); ++i)
    EXPECT_FLOAT_EQ(a[i], b[i]);
}

#endif
