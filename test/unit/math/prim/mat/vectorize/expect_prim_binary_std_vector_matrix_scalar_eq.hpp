#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_PRIM_BINARY_STD_VECTOR_MATRIX_SCALAR_EQ
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_PRIM_BINARY_STD_VECTOR_MATRIX_SCALAR_EQ

#include <stan/math/prim/scal/fun/is_nan.hpp>
#include <Eigen/Dense>
#include <vector>
#include <gtest/gtest.h>

template <typename F, typename input_t, typename matrix_t>
void expect_prim_binary_std_vector_matrix_scalar_eq( 
const std::vector<matrix_t>& input_mv, input_t input) {
  using stan::math::is_nan;
  using std::vector;

  vector<matrix_t> fd = 
  F::template apply<vector<matrix_t> >(input_mv, input);
  EXPECT_EQ(input_mv.size(), fd.size());
  for (size_t i = 0; i < input_mv.size(); ++i) {
    EXPECT_EQ(input_mv[i].size(), fd[i].size());
    for (int j = 0; j < input_mv[i].size(); ++j) {
      double exp_v = F::apply_base(input_mv[i](j), input);
      if (is_nan(exp_v) && is_nan(fd[i](j))) continue;
      EXPECT_FLOAT_EQ(exp_v, fd[i](j));
    }
  }    
}
#endif
