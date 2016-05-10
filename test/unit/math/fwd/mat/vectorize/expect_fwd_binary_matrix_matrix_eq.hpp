#ifndef TEST_UNIT_MATH_FWD_MAT_VECTORIZE_EXPECT_FWD_BINARY_MATRIX_MATRIX_EQ
#define TEST_UNIT_MATH_FWD_MAT_VECTORIZE_EXPECT_FWD_BINARY_MATRIX_MATRIX_EQ

#include <stan/math/fwd/core/fvar.hpp>
#include <test/unit/math/fwd/mat/vectorize/build_fwd_binary_matrix.hpp>
#include <test/unit/math/fwd/mat/vectorize/expect_val_deriv_eq.hpp>
#include <Eigen/Dense>

template <typename F, typename FV, typename matrix_t1, typename matrix_t2>
void expect_fwd_binary_matrix_matrix_eq(const matrix_t1& template_m1, 
const matrix_t2& template_m2, bool seed_one = 1, bool seed_two = 1) {

  typedef Eigen::Matrix<FV, matrix_t1::RowsAtCompileTime, 
  matrix_t1::ColsAtCompileTime> result_matrix_t;

  for (int i = 0; i < template_m1.size(); ++i) {
    matrix_t1 input_m1;
    if (seed_one)
      input_m1 = build_fwd_binary_matrix1<F>(template_m1, i);
    else
      input_m1 = build_fwd_binary_matrix1<F>(template_m1);
    matrix_t2 input_m2;
    if (seed_two)
      input_m2 = build_fwd_binary_matrix2<F>(template_m2, i);
    else
      input_m2 = build_fwd_binary_matrix2<F>(template_m2);
    result_matrix_t fa = F::template apply<result_matrix_t>(
    input_m1, input_m2);
    EXPECT_EQ(input_m1.size(), fa.size());
    EXPECT_EQ(input_m2.size(), fa.size());
    expect_val_deriv_eq(F::apply_base(input_m1(i), input_m2(i)),
    fa(i));
  }   

  if (template_m1.rows() > 1 && template_m1.cols() > 1) {
    matrix_t1 input_m1 = build_fwd_binary_matrix1<F>(template_m1);
    matrix_t2 input_m2 = build_fwd_binary_matrix2<F>(template_m2);
 
    result_matrix_t fb = F::template apply<result_matrix_t>(
    input_m1.block(1, 1, 1, 1), input_m2.block(1, 1, 1, 1)); 
    expect_val_deriv_eq(F::apply_base(input_m1(1, 1), 
    input_m2(1, 1)), fb(0, 0));  
  }
}
#endif
