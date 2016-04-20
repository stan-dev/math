#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_PRIM_BINARY_STD_VECTOR_MATRIX_STD_VECTOR_MATRIX_EQ
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_PRIM_BINARY_STD_VECTOR_MATRIX_STD_VECTOR_MATRIX_EQ

#include <stan/math/prim/scal/fun/is_nan.hpp>
#include <Eigen/Dense>
#include <vector>
#include <gtest/gtest.h>


template <typename F, typename matrix_t>
void expect_prim_binary_std_vector_matrix_std_vector_matrix_eq(const
std::vector<matrix_t>& input_mv, int m_size) {
  using stan::math::is_nan;
  using std::vector;

  std::vector<matrix_t> fa = F::template apply<vector<matrix_t> >(
  input_mv, input_mv);
  EXPECT_EQ(input_mv.size(), fa.size());
  for (size_t i = 0; i < fa.size(); ++i) {
    EXPECT_EQ(m_size, fa[i].size());
    for (int j = 0; j < fa[i].size(); ++j) {
      double exp_v = F::apply_base(input_mv[i](j), input_mv[i](j));
      if (is_nan(exp_v) && is_nan(fa[i](j))) continue;
      EXPECT_FLOAT_EQ(exp_v, fa[i](j));
    }
  }    
}
#endif
