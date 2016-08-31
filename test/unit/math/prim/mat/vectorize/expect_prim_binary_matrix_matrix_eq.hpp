#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_PRIM_BINARY_MATRIX_MATRIX_EQ
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_PRIM_BINARY_MATRIX_MATRIX_EQ

#include <test/unit/math/prim/mat/vectorize/expect_val_eq.hpp>
#include <Eigen/Dense>
#include <gtest/gtest.h>


template <typename F, typename matrix_t>
void expect_prim_binary_matrix_matrix_eq(const matrix_t& input_m1, 
const matrix_t& input_m2) {

  matrix_t fa = F::template apply<matrix_t>(input_m1, input_m2);
  EXPECT_EQ(input_m1.size(), fa.size());
  EXPECT_EQ(input_m2.size(), fa.size());
  for (int i = 0; i < input_m1.size(); ++i) {
    double exp_v = F::apply_base(input_m1(i), input_m2(i));
    expect_val_eq(exp_v, fa(i));
  }    
}
#endif
