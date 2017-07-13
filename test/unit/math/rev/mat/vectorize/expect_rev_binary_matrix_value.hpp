#ifndef TEST_UNIT_MATH_REV_MAT_VECTORIZE_EXPECT_REV_BINARY_MATRIX_VALUE
#define TEST_UNIT_MATH_REV_MAT_VECTORIZE_EXPECT_REV_BINARY_MATRIX_VALUE

#include <stan/math/rev/core/var.hpp>
#include <test/unit/math/rev/mat/vectorize/build_binary_vector.hpp>
#include <test/unit/math/rev/mat/vectorize/build_rev_binary_matrix.hpp>
#include <test/unit/math/prim/mat/vectorize/build_template_matrix.hpp>
#include <test/unit/math/rev/mat/vectorize/expect_binary_val_deriv_eq.hpp>
#include <vector>
#include <Eigen/Dense>

template <typename F, typename input_t1, typename input_t2, int R, int C>
void expect_rev_binary_scalar_matrix_eq(
  const std::vector<input_t1>& template_scalar_v,
  const Eigen::Matrix<input_t2, R, C>& template_m) {

  using std::vector;
  using stan::math::var;

  typedef Eigen::Matrix<input_t2, R, C> input_matrix_t;
  typedef Eigen::Matrix<var, R, C> result_matrix_t;

  for (size_t i = 0; i < template_scalar_v.size(); ++i) {
    for (int j = 0; j < template_m.size(); ++j) {
      vector<input_t1> input_va1 = build_binary_vector1<F>(
        template_scalar_v);
      vector<input_t1> input_va2 = build_binary_vector1<F>(
        template_scalar_v);
      vector<input_t2> input_vb1 = build_binary_vector2<F>(
        vector<input_t2>());
      vector<input_t2> input_vb2 = build_binary_vector2<F>(
        vector<input_t2>());
      input_matrix_t input_m1 = build_rev_binary_matrix(input_vb1[i],
                                                        template_m);
      input_matrix_t input_m2 = build_rev_binary_matrix(input_vb2[i],
                                                        template_m);
      result_matrix_t fa = F::template apply<result_matrix_t>(
        input_va2[i], input_m2);
      EXPECT_EQ(input_m2.size(), fa.size());
      expect_binary_val_deriv_eq(F::apply_base(input_va1[i], input_m1(j)),
                                 input_va1[i], input_m1(j), fa(j),
                                 input_va2[i], input_m2(j));
    }
  }
}

template <typename F, typename input_t1, typename input_t2, int R, int C>
void expect_rev_binary_matrix_scalar_eq(
  const Eigen::Matrix<input_t1, R, C>& template_m,
  const std::vector<input_t2>& template_scalar_v) {

  using std::vector;
  using stan::math::var;

  typedef Eigen::Matrix<input_t1, R, C> input_matrix_t;
  typedef Eigen::Matrix<var, R, C> result_matrix_t;

  for (size_t i = 0; i < template_scalar_v.size(); ++i) {
    for (int j = 0; j < template_m.size(); ++j) {
      vector<input_t1> input_va1 = build_binary_vector1<F>(
        vector<input_t1>());
      vector<input_t1> input_va2 = build_binary_vector1<F>(
        vector<input_t1>());
      vector<input_t2> input_vb1 = build_binary_vector2<F>(
        template_scalar_v);
      vector<input_t2> input_vb2 = build_binary_vector2<F>(
        template_scalar_v);
      input_matrix_t input_m1 = build_rev_binary_matrix(input_va1[i],
        template_m);
      input_matrix_t input_m2 = build_rev_binary_matrix(input_va2[i],
        template_m);
      result_matrix_t fa = F::template apply<result_matrix_t>(
        input_m2, input_vb2[i]);
      EXPECT_EQ(input_m2.size(), fa.size());
      expect_binary_val_deriv_eq(F::apply_base(input_m1(j), input_vb1[i]),
                                 input_m1(j), input_vb1[i], fa(j),
                                 input_m2(j), input_vb2[i]);
    }
  }
}

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

template <typename F, typename input_t1, typename input_t2, int R, int C>
void expect_rev_binary_scalar_std_vector_matrix_eq(
const std::vector<input_t1>& template_scalar_v,
const Eigen::Matrix<input_t2, R, C>& template_m) {
  using std::vector;
  using stan::math::var;

  typedef Eigen::Matrix<input_t2, R, C> input_matrix_t;
  typedef Eigen::Matrix<var, R, C> result_matrix_t;

  for (size_t i = 0; i < template_scalar_v.size(); ++i) {
    for (size_t j = 0; j < 2; ++j) {
      for (int k = 0; k < template_m.size(); ++k) {
        vector<input_t1> input_va1 = build_binary_vector1<F>(
          template_scalar_v);
        vector<input_t1> input_va2 = build_binary_vector1<F>(
          template_scalar_v);
        vector<input_t2> input_vb1 = build_binary_vector2<F>(
          vector<input_t2>());
        vector<input_t2> input_vb2 = build_binary_vector2<F>(
          vector<input_t2>());
        vector<input_matrix_t> input_mv1;
        input_mv1.push_back(build_rev_binary_matrix(input_vb1[i],
                                                    template_m));
        input_mv1.push_back(build_rev_binary_matrix(input_vb1[i],
                                                    template_m));
        vector<input_matrix_t> input_mv2;
        input_mv2.push_back(build_rev_binary_matrix(input_vb2[i],
                                                    template_m));
        input_mv2.push_back(build_rev_binary_matrix(input_vb2[i],
                                                    template_m));
        vector<result_matrix_t> fa = F::template
          apply<vector<result_matrix_t> >(input_va2[i], input_mv2);
        EXPECT_EQ(input_mv2.size(), fa.size());
        EXPECT_EQ(input_mv2[j].size(), fa[j].size());
        expect_binary_val_deriv_eq(F::apply_base(input_va1[i],
                                                 input_mv1[j](k)),
                                   input_va1[i], input_mv1[j](k),
                                   fa[j](k), input_va2[i], input_mv2[j](k));
      }
    }
  }
}

template <typename F, typename input_t1, typename input_t2, int R, int C>
void expect_rev_binary_std_vector_matrix_scalar_eq(
const Eigen::Matrix<input_t1, R, C>& template_m,
const std::vector<input_t2>& template_scalar_v) {
  using std::vector;
  using stan::math::var;

  typedef Eigen::Matrix<input_t1, R, C> input_matrix_t;
  typedef Eigen::Matrix<var, R, C> result_matrix_t;

  for (size_t i = 0; i < template_scalar_v.size(); ++i) {
    for (size_t j = 0; j < 2; ++j) {
      for (int k = 0; k < template_m.size(); ++k) {
        vector<input_t1> input_va1 = build_binary_vector1<F>(
          vector<input_t1>());
        vector<input_t1> input_va2 = build_binary_vector1<F>(
          vector<input_t1>());
        vector<input_t2> input_vb1 = build_binary_vector2<F>(
          template_scalar_v);
        vector<input_t2> input_vb2 = build_binary_vector2<F>(
          template_scalar_v);
        vector<input_matrix_t> input_mv1;
        input_mv1.push_back(build_rev_binary_matrix(input_va1[i],
                                                    template_m));
        input_mv1.push_back(build_rev_binary_matrix(input_va1[i],
                                                    template_m));
        vector<input_matrix_t> input_mv2;
        input_mv2.push_back(build_rev_binary_matrix(input_va2[i],
                                                    template_m));
        input_mv2.push_back(build_rev_binary_matrix(input_va2[i],
                                                    template_m));
        vector<result_matrix_t> fa = F::template
          apply<vector<result_matrix_t> >(input_mv2, input_vb2[i]);
        EXPECT_EQ(input_mv2.size(), fa.size());
        EXPECT_EQ(input_mv2[j].size(), fa[j].size());
        expect_binary_val_deriv_eq(F::apply_base(input_mv1[j](k),
                                                 input_vb1[i]),
                                   input_mv1[j](k), input_vb1[i],
                                   fa[j](k), input_mv2[j](k), input_vb2[i]);
      }
    }
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

template <typename F, int R, int C>
void expect_rev_binary_scalar_matrix_all_eq(
  const std::vector<int>& int_template_v,
  const std::vector<double>& d_template_v,
  const std::vector<stan::math::var>& var_template_v,
  const Eigen::Matrix<double, R, C> d_scalar_template_m,
  const Eigen::Matrix<stan::math::var, R, C> var_scalar_template_m) {

  expect_rev_binary_scalar_matrix_eq<F>(int_template_v,
                                        var_scalar_template_m);
  expect_rev_binary_matrix_scalar_eq<F>(var_scalar_template_m,
                                        int_template_v);

  expect_rev_binary_scalar_matrix_eq<F>(var_template_v,
                                        d_scalar_template_m);
  expect_rev_binary_matrix_scalar_eq<F>(d_scalar_template_m,
                                        var_template_v);
  expect_rev_binary_scalar_matrix_eq<F>(d_template_v,
                                        var_scalar_template_m);
  expect_rev_binary_matrix_scalar_eq<F>(var_scalar_template_m,
                                        d_template_v);

  expect_rev_binary_scalar_matrix_eq<F>(var_template_v,
                                        var_scalar_template_m);
  expect_rev_binary_matrix_scalar_eq<F>(var_scalar_template_m,
                                        var_template_v);
}

template <typename F, int R, int C>
void expect_rev_binary_scalar_std_vector_matrix_all_eq(
  const std::vector<int>& int_template_v,
  const std::vector<double>& d_template_v,
  const std::vector<stan::math::var>& var_template_v,
  const Eigen::Matrix<double, R, C> d_scalar_template_m,
  const Eigen::Matrix<stan::math::var, R, C> var_scalar_template_m) {

  expect_rev_binary_scalar_std_vector_matrix_eq<F>(int_template_v,
                                                   var_scalar_template_m);
  expect_rev_binary_std_vector_matrix_scalar_eq<F>(var_scalar_template_m,
                                                   int_template_v);

  expect_rev_binary_scalar_std_vector_matrix_eq<F>(var_template_v,
                                                   d_scalar_template_m);
  expect_rev_binary_std_vector_matrix_scalar_eq<F>(d_scalar_template_m,
                                                   var_template_v);
  expect_rev_binary_scalar_std_vector_matrix_eq<F>(d_template_v,
                                                   var_scalar_template_m);
  expect_rev_binary_std_vector_matrix_scalar_eq<F>(var_scalar_template_m,
                                                   d_template_v);

  expect_rev_binary_scalar_std_vector_matrix_eq<F>(var_template_v,
                                                   var_scalar_template_m);
  expect_rev_binary_std_vector_matrix_scalar_eq<F>(var_scalar_template_m,
                                                   var_template_v);
}

template <typename F, typename matrix_t, int R, int C>
void expect_rev_binary_matrix_value(
  Eigen::Matrix<matrix_t, R, C> template_m) {

  using stan::math::var;
  using std::vector;

  vector<int> int_template_v;
  vector<double> d_template_v;
  vector<var> var_template_v;

  typedef Eigen::Matrix<double, R, C> d_matrix_t;
  typedef Eigen::Matrix<var, R, C> var_matrix_t;

  d_matrix_t d_scalar_template_m;
  d_matrix_t d_template_m;
  var_matrix_t var_scalar_template_m;
  var_matrix_t var_template_m;

  d_scalar_template_m = build_template_matrix(d_scalar_template_m, 3, 5);
  d_template_m = build_template_matrix(d_scalar_template_m, 
                                       3, F::valid_inputs1().size());
  var_scalar_template_m = build_template_matrix(var_scalar_template_m, 
                                                3, 5);
  var_template_m = build_template_matrix(var_scalar_template_m, 
                                         3, F::valid_inputs1().size());
  expect_rev_binary_scalar_matrix_all_eq<F>(int_template_v, d_template_v,
                                            var_template_v,
                                            d_scalar_template_m,
                                            var_scalar_template_m);

  expect_rev_binary_matrix_matrix_eq<F>(var_template_m, d_template_m);
  expect_rev_binary_matrix_matrix_eq<F>(d_template_m, var_template_m);
  expect_rev_binary_matrix_matrix_eq<F>(var_template_m, var_template_m);

  expect_rev_binary_scalar_std_vector_matrix_all_eq<F>(
    int_template_v, d_template_v, var_template_v,
    d_scalar_template_m, var_scalar_template_m);

  expect_rev_binary_std_vector_matrix_std_vector_matrix_eq<F>(
    var_template_m, d_template_m);
  expect_rev_binary_std_vector_matrix_std_vector_matrix_eq<F>(
    d_template_m, var_template_m);
  expect_rev_binary_std_vector_matrix_std_vector_matrix_eq<F>(
    var_template_m, var_template_m);
}

#endif
