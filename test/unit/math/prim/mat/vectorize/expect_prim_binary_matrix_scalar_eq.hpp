#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_PRIM_BINARY_MATRIX_SCALAR_EQ
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_PRIM_BINARY_MATRIX_SCALAR_EQ

#include <test/unit/math/prim/mat/vectorize/expect_val_eq.hpp>
#include <Eigen/Dense>
#include <gtest/gtest.h>

template <typename F, typename input_t, typename matrix_t>
void expect_prim_binary_matrix_scalar_eq(const matrix_t& input_m,
                                         input_t input) {
  matrix_t fa = F::template apply<matrix_t>(input_m, input);
  EXPECT_EQ(input_m.size(), fa.size());
  for (int i = 0; i < fa.size(); ++i) {
    double exp_v = F::apply_base(input_m(i), input);
    expect_val_eq(exp_v, fa(i));
  }
}
#endif
