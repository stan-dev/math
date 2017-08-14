#ifndef TEST_UNIT_MATH_REV_MAT_VECTORIZE_EXPECT_REV_BINARY_MATRIX_VALUE
#define TEST_UNIT_MATH_REV_MAT_VECTORIZE_EXPECT_REV_BINARY_MATRIX_VALUE

#include <stan/math/rev/core/var.hpp>
#include <test/unit/math/rev/mat/vectorize/build_binary_vector.hpp>
#include <test/unit/math/rev/mat/vectorize/build_rev_binary_matrix.hpp>
#include <test/unit/math/prim/mat/vectorize/build_template_matrix.hpp>
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
                               input_ma1(i), input_mb1(i), fa(i),
                               input_ma2(i), input_mb2(i));
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

template <typename F, typename matrix_t1, typename matrix_t2>
void expect_rev_binary_std_vector_matrix_std_vector_matrix_eq(
const matrix_t1& template_m1, const matrix_t2& template_m2) {
  using std::vector;
  using stan::math::var;

  typedef Eigen::Matrix<var, matrix_t1::RowsAtCompileTime,
  matrix_t1::ColsAtCompileTime> result_matrix_t;

  for (size_t i = 0; i < 2; ++i) {
    for (int j = 0; j < template_m1.size(); ++j) {
      vector<matrix_t1> input_mva1;
      input_mva1.push_back(build_rev_binary_matrix1<F>(template_m1));
      input_mva1.push_back(build_rev_binary_matrix1<F>(template_m1));
      vector<matrix_t1> input_mva2;
      input_mva2.push_back(build_rev_binary_matrix1<F>(template_m1));
      input_mva2.push_back(build_rev_binary_matrix1<F>(template_m1));
      vector<matrix_t2> input_mvb1;
      input_mvb1.push_back(build_rev_binary_matrix2<F>(template_m2));
      input_mvb1.push_back(build_rev_binary_matrix2<F>(template_m2));
      vector<matrix_t2> input_mvb2;
      input_mvb2.push_back(build_rev_binary_matrix2<F>(template_m2));
      input_mvb2.push_back(build_rev_binary_matrix2<F>(template_m2));
      vector<result_matrix_t> fa = F::template
      apply<vector<result_matrix_t> >(input_mva2, input_mvb2);
      EXPECT_EQ(input_mva2.size(), fa.size());
      EXPECT_EQ(input_mvb2.size(), fa.size());
      EXPECT_EQ(input_mva2[i].size(), fa[i].size());
      EXPECT_EQ(input_mvb2[i].size(), fa[i].size());
      expect_binary_val_deriv_eq(F::apply_base(input_mva1[i](j),
                                               input_mvb1[i](j)),
                                 input_mva1[i](j), input_mvb1[i](j),
                                 fa[i](j), input_mva2[i](j),
                                 input_mvb2[i](j));
    }
  }
}

template <typename F, int R, int C> void
expect_rev_binary_matrix_value(Eigen::Matrix<double, R, C> template_m) {

  using stan::math::var;
  using std::vector;

  vector<int> int_template_v;
  vector<double> d_template_v;
  vector<var> var_template_v;

  typedef Eigen::Matrix<double, R, C> d_matrix_t;
  typedef Eigen::Matrix<var, R, C> var_matrix_t;

  d_matrix_t d_template_m;
  var_matrix_t var_template_m;

  d_template_m = build_template_matrix(d_template_m,
                                       F::valid_inputs1().size(), 3);
  var_template_m = build_template_matrix(var_template_m,
                                         F::valid_inputs1().size(), 3);

  expect_rev_binary_matrix_matrix_eq<F>(var_template_m, d_template_m);
  expect_rev_binary_matrix_matrix_eq<F>(d_template_m, var_template_m);
  expect_rev_binary_matrix_matrix_eq<F>(var_template_m, var_template_m);

  expect_rev_binary_std_vector_matrix_std_vector_matrix_eq<F>(
    var_template_m, d_template_m);
  expect_rev_binary_std_vector_matrix_std_vector_matrix_eq<F>(
    d_template_m, var_template_m);
  expect_rev_binary_std_vector_matrix_std_vector_matrix_eq<F>(
    var_template_m, var_template_m);
}

#endif
