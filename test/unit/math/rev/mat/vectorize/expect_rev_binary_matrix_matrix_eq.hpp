#ifndef TEST_UNIT_MATH_REV_MAT_VECTORIZE_EXPECT_REV_BINARY_MATRIX_MATRIX_EQ
#define TEST_UNIT_MATH_REV_MAT_VECTORIZE_EXPECT_REV_BINARY_MATRIX_MATRIX_EQ

#include <stan/math/rev/core/var.hpp>
#include <test/unit/math/rev/mat/vectorize/build_rev_binary_matrix.hpp>
#include <test/unit/math/rev/mat/vectorize/expect_binary_val_deriv_eq.hpp>
#include <vector>
#include <Eigen/Dense>

template <typename F, typename matrix_t1, typename matrix_t2>
void expect_rev_binary_matrix_matrix_eq(const matrix_t1& template_m1, 
const matrix_t2& template_m2) {
  using std::vector;
  using stan::math::var;

  typedef Eigen::Matrix<var, matrix_t1::RowsAtCompileTime, 
  matrix_t1::ColsAtCompileTime> result_matrix_t;

  for (int i = 0; i < template_m1.size(); ++i) {
    matrix_t1 input_ma1 = build_rev_binary_matrix1<F>(template_m1);
    matrix_t1 input_ma2 = build_rev_binary_matrix1<F>(template_m1);
    matrix_t2 input_mb1 = build_rev_binary_matrix2<F>(template_m2);
    matrix_t2 input_mb2 = build_rev_binary_matrix2<F>(template_m2);
    result_matrix_t fa = F::template apply<result_matrix_t>(
    input_ma2, input_mb2);
    EXPECT_EQ(input_ma2.size(), fa.size());
    EXPECT_EQ(input_mb2.size(), fa.size());
    expect_binary_val_deriv_eq(F::apply_base(input_ma1(i), input_mb1(i)),
    input_ma1(i), input_mb1(i), fa(i), input_ma2(i), input_mb2(i));
  }   

  if (template_m1.rows() > 1 && template_m1.cols() > 1) {
    matrix_t1 input_ma1 = build_rev_binary_matrix1<F>(template_m1);
    matrix_t1 input_ma2 = build_rev_binary_matrix1<F>(template_m1);
    matrix_t2 input_mb1 = build_rev_binary_matrix2<F>(template_m2);
    matrix_t2 input_mb2 = build_rev_binary_matrix2<F>(template_m2);
 
    result_matrix_t fb = F::template apply<result_matrix_t>(
    input_ma2.block(1, 1, 1, 1), input_mb2.block(1, 1, 1, 1)); 
    expect_binary_val_deriv_eq(F::apply_base(input_ma1(1, 1), 
    input_mb1(1, 1)), input_ma1(1, 1), input_mb1(1, 1), 
    fb(0, 0), input_ma2(1, 1), input_mb2(1, 1));  
  }
}
#endif
