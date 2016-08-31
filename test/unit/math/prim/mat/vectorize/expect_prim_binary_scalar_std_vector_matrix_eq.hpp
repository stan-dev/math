#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_PRIM_BINARY_SCALAR_STD_VECTOR_MATRIX_EQ
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_PRIM_BINARY_SCALAR_STD_VECTOR_MATRIX_EQ

#include <test/unit/math/prim/mat/vectorize/expect_val_eq.hpp>
#include <Eigen/Dense>
#include <vector>
#include <gtest/gtest.h>

template <typename F, typename input_t, typename matrix_t>
void expect_prim_binary_scalar_std_vector_matrix_eq(input_t input, 
const std::vector<matrix_t>& input_mv) {
  using std::vector;

  vector<matrix_t> fd=F::template apply<vector<matrix_t> >(input,input_mv);
  EXPECT_EQ(input_mv.size(), fd.size());
  for (size_t i = 0; i < input_mv.size(); ++i) {
    EXPECT_EQ(input_mv[i].size(), fd[i].size());
    for (int j = 0; j < input_mv[i].size(); ++j) {
      double exp_v = F::apply_base(input, input_mv[i](j));
      expect_val_eq(exp_v, fd[i](j));
    }
  }    
}
#endif
