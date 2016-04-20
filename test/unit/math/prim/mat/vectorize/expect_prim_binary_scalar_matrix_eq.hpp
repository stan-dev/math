#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_PRIM_BINARY_SCALAR_MATRIX_EQ
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_PRIM_BINARY_SCALAR_MATRIX_EQ

#include <stan/math/prim/scal/fun/is_nan.hpp>
#include <Eigen/Dense>
#include <gtest/gtest.h>


template <typename F, typename input_t, typename matrix_t>
void expect_prim_binary_scalar_matrix_eq(input_t input, const matrix_t&
input_m) {
  using stan::math::is_nan;

  matrix_t fa = F::template apply<matrix_t>(input, input_m);
  EXPECT_EQ(input_m.size(), fa.size());
  for (int i = 0; i < fa.size(); ++i) {
    double exp_v = F::apply_base(input, input_m(i));
    if (is_nan(exp_v) && is_nan(fa(i))) continue;
    EXPECT_FLOAT_EQ(exp_v, fa(i));
  }    

  fa = F::template apply<matrix_t>(input_m, input);
  EXPECT_EQ(input_m.size(), fa.size());
  for (int i = 0; i < fa.size(); ++i) {
    double exp_v = F::apply_base(input_m(i), input);
    if (is_nan(exp_v) && is_nan(fa(i))) continue;
    EXPECT_FLOAT_EQ(exp_v, fa(i));
  }    
}
#endif
