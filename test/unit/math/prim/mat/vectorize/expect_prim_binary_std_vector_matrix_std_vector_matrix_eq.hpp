#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_PRIM_BINARY_STD_VECTOR_MATRIX_STD_VECTOR_MATRIX_EQ
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_PRIM_BINARY_STD_VECTOR_MATRIX_STD_VECTOR_MATRIX_EQ

#include <test/unit/math/prim/mat/vectorize/expect_val_eq.hpp>
#include <Eigen/Dense>
#include <vector>
#include <gtest/gtest.h>


template <typename F, typename matrix_t> void
expect_prim_binary_std_vector_matrix_std_vector_matrix_eq(
  const std::vector<matrix_t>& input_mv1,
  const std::vector<matrix_t>& input_mv2) {

  using std::vector;

  std::vector<matrix_t> fa = F::template apply<vector<matrix_t> >(
    input_mv1, input_mv2);
  EXPECT_EQ(input_mv1.size(), fa.size());
  EXPECT_EQ(input_mv2.size(), fa.size());
  for (size_t i = 0; i < input_mv1.size(); ++i) {
    EXPECT_EQ(input_mv1[i].size(), fa[i].size());
    EXPECT_EQ(input_mv2[i].size(), fa[i].size());
    for (int j = 0; j < input_mv1[i].size(); ++j) {
      double exp_v = F::apply_base(input_mv1[i](j), input_mv2[i](j));
      expect_val_eq(exp_v, fa[i](j));
    }
  }
}
#endif
