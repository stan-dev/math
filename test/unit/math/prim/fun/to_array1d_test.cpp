#include <stan/math/prim.hpp>
#include <test/unit/math/prim/fun/expect_matrix_eq.hpp>
#include <gtest/gtest.h>
#include <vector>
#include <stdexcept>

using stan::math::to_array_1d;

TEST(MathMatrix, to_array_1d_matrix) {
  Eigen::MatrixXd a(3, 3);
  a << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  std::vector<double> a_correct{1, 4, 7, 2, 5, 8, 3, 6, 9};
  std::vector<double> a_res = to_array_1d(a);
  expect_std_vector_eq(a_res, a_correct);
}

TEST(MathMatrix, to_array_1d_matrix_block) {
  Eigen::MatrixXd a(3, 3);
  a << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  std::vector<double> a_correct{2, 5, 3, 6};
  std::vector<double> a_res = to_array_1d(a.block(0, 1, 2, 2));
  expect_std_vector_eq(a_res, a_correct);
}
