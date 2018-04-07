#ifndef TEST_UNIT_MATH_FWD_MAT_VECTORIZE_EXPECT_MIX_BINARY_MATRIX_VALUE
#define TEST_UNIT_MATH_FWD_MAT_VECTORIZE_EXPECT_MIX_BINARY_MATRIX_VALUE

#include <stan/math/rev/core/var.hpp>
#include <stan/math/fwd/core/fvar.hpp>
#include <test/unit/math/rev/mat/vectorize/build_binary_vector.hpp>
#include <test/unit/math/fwd/mat/vectorize/build_fwd_binary_matrix.hpp>
#include <test/unit/math/prim/mat/vectorize/build_template_matrix.hpp>
#include <test/unit/math/mix/mat/vectorize/expect_binary_val_deriv_eq.hpp>
#include <Eigen/Dense>
#include <vector>

template <typename F, typename FV, typename matrix_t1, typename matrix_t2>
void expect_mix_binary_matrix_matrix_eq(const matrix_t1& template_m1,
                                        const matrix_t2& template_m2,
                                        bool seed_one = 1,
                                        bool seed_two = 1) {
  typedef Eigen::Matrix<FV, matrix_t1::RowsAtCompileTime,
                        matrix_t1::ColsAtCompileTime> result_matrix_t;

  for (int i = 0; i < template_m1.size(); ++i) {
    matrix_t1 input_ma1;
    matrix_t1 input_ma2;
    if (seed_one) {
      input_ma1 = build_fwd_binary_matrix1<F>(template_m1, i);
      input_ma2 = build_fwd_binary_matrix1<F>(template_m1, i);
    } else {
      input_ma1 = build_fwd_binary_matrix1<F>(template_m1);
      input_ma2 = build_fwd_binary_matrix1<F>(template_m1);
    }
    matrix_t2 input_mb1;
    matrix_t2 input_mb2;
    if (seed_two) {
      input_mb1 = build_fwd_binary_matrix2<F>(template_m2, i);
      input_mb2 = build_fwd_binary_matrix2<F>(template_m2, i);
    } else {
      input_mb1 = build_fwd_binary_matrix2<F>(template_m2);
      input_mb2 = build_fwd_binary_matrix2<F>(template_m2);
    }
    result_matrix_t fa = F::template apply<result_matrix_t>(input_ma2,
                                                            input_mb2);
    EXPECT_EQ(input_ma2.size(), fa.size());
    EXPECT_EQ(input_mb2.size(), fa.size());
    expect_binary_val_deriv_eq(F::apply_base(input_ma1(i), input_mb1(i)),
                               input_ma1(i), input_mb1(i),
                               fa(i), input_ma2(i), input_mb2(i));
  }

  if (template_m1.rows() > 1 && template_m1.cols() > 1) {
    matrix_t1 input_ma1 = build_fwd_binary_matrix1<F>(template_m1);
    matrix_t1 input_ma2 = build_fwd_binary_matrix1<F>(template_m1);
    matrix_t2 input_mb1 = build_fwd_binary_matrix2<F>(template_m2);
    matrix_t2 input_mb2 = build_fwd_binary_matrix2<F>(template_m2);

    result_matrix_t fb = F::template apply<result_matrix_t>(
      input_ma2.block(1, 1, 1, 1), input_mb2.block(1, 1, 1, 1));
    expect_binary_val_deriv_eq(F::apply_base(input_ma1(1, 1),
                                             input_mb1(1, 1)),
                               input_ma1(1, 1), input_mb1(1, 1),
                               fb(0, 0), input_ma2(1, 1), input_mb2(1, 1));
  }
}

template <typename F, typename FV, typename matrix_t1, typename matrix_t2>
void expect_mix_binary_std_vector_matrix_std_vector_matrix_eq(
    const matrix_t1& template_m1, const matrix_t2& template_m2,
    bool seed_one = 1, bool seed_two = 1) {
  using std::vector;
  typedef Eigen::Matrix<FV, matrix_t1::RowsAtCompileTime,
                        matrix_t1::ColsAtCompileTime> result_matrix_t;

  const size_t num_v = 2;
  for (size_t i = 0; i < num_v; ++i) {
    for (int j = 0; j < template_m1.size(); ++j) {
      vector<matrix_t1> input_mva1;
      vector<matrix_t1> input_mva2;
      for (size_t k = 0; k < num_v; ++k) {
        if (seed_one && k == i) {
          input_mva1.push_back(build_fwd_binary_matrix1<F>(template_m1, j));
          input_mva2.push_back(build_fwd_binary_matrix1<F>(template_m1, j));
        } else {
          input_mva1.push_back(build_fwd_binary_matrix1<F>(template_m1));
          input_mva2.push_back(build_fwd_binary_matrix1<F>(template_m1));
        }
      }
      vector<matrix_t2> input_mvb1;
      vector<matrix_t2> input_mvb2;
      for (size_t k = 0; k < num_v; ++k) {
        if (seed_one && k == i) {
          input_mvb1.push_back(build_fwd_binary_matrix2<F>(template_m2, j));
          input_mvb2.push_back(build_fwd_binary_matrix2<F>(template_m2, j));
        } else {
          input_mvb1.push_back(build_fwd_binary_matrix2<F>(template_m2));
          input_mvb2.push_back(build_fwd_binary_matrix2<F>(template_m2));
        }
      }
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

template <typename F, typename FV, int R, int C>
void expect_mix_binary_matrix_value(
    Eigen::Matrix<double, R, C> template_m) {
  using std::vector;

  vector<int> int_template_v;
  vector<double> d_template_v;
  vector<FV> mix_template_v;

  typedef Eigen::Matrix<double, R, C> d_matrix_t;
  typedef Eigen::Matrix<FV, R, C> mix_matrix_t;

  d_matrix_t d_template_m;
  mix_matrix_t mix_template_m;

  d_template_m = build_template_matrix(d_template_m,
                                       F::valid_inputs1().size(), 3);
  mix_template_m = build_template_matrix(mix_template_m,
                                         F::valid_inputs1().size(), 3);

  expect_mix_binary_matrix_matrix_eq<F, FV>(mix_template_m, d_template_m);
  expect_mix_binary_matrix_matrix_eq<F, FV>(d_template_m, mix_template_m);
  expect_mix_binary_matrix_matrix_eq<F, FV>(mix_template_m,
                                            mix_template_m, 1, 0);
  expect_mix_binary_matrix_matrix_eq<F, FV>(mix_template_m,
                                            mix_template_m, 0, 1);
  expect_mix_binary_matrix_matrix_eq<F, FV>(mix_template_m,
                                            mix_template_m, 1, 1);

  expect_mix_binary_std_vector_matrix_std_vector_matrix_eq<F, FV>(
      mix_template_m, d_template_m);
  expect_mix_binary_std_vector_matrix_std_vector_matrix_eq<F, FV>(
      d_template_m, mix_template_m);
  expect_mix_binary_std_vector_matrix_std_vector_matrix_eq<F, FV>(
      mix_template_m, mix_template_m, 1, 0);
  expect_mix_binary_std_vector_matrix_std_vector_matrix_eq<F, FV>(
      mix_template_m, mix_template_m, 0, 1);
  expect_mix_binary_std_vector_matrix_std_vector_matrix_eq<F, FV>(
      mix_template_m, mix_template_m, 1, 1);
}
#endif
