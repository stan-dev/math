#ifndef TEST_UNIT_MATH_FWD_MAT_VECTORIZE_EXPECT_FWD_BINARY_MATRIX_VALUE
#define TEST_UNIT_MATH_FWD_MAT_VECTORIZE_EXPECT_FWD_BINARY_MATRIX_VALUE

#include <stan/math/fwd/core/fvar.hpp>
#include <test/unit/math/fwd/mat/vectorize/build_binary_vector.hpp>
#include <test/unit/math/fwd/mat/vectorize/build_fwd_binary_matrix.hpp>
#include <test/unit/math/prim/mat/vectorize/build_template_matrix.hpp>
#include <test/unit/math/fwd/mat/vectorize/expect_val_deriv_eq.hpp>
#include <vector>
#include <Eigen/Dense>

template <typename F, typename FV, typename input_t1, typename input_t2, 
          int R, int C> 
void expect_fwd_binary_scalar_matrix_eq(
  const std::vector<input_t1>& template_scalar_v, 
  const Eigen::Matrix<input_t2, R, C>& template_m, 
  bool seed_one = 1, bool seed_two = 1) {

  using std::vector;

  typedef Eigen::Matrix<input_t2, R, C> input_matrix_t;
  typedef Eigen::Matrix<FV, R, C> result_matrix_t;

  for (size_t i = 0; i < template_scalar_v.size(); ++i) {
    input_t1 input_1 = 0;
    if (seed_one)
      input_1 = build_binary_vector1<F>(template_scalar_v, i)[i];
    else
      input_1 = build_binary_vector1<F>(template_scalar_v)[i];
    for (int j = 0; j < template_m.size(); ++j) {
      input_matrix_t input_m(template_m.rows(), template_m.cols()); 
      if (seed_two)
        input_m = build_fwd_binary_matrix2<F>(i, template_m, j);
      else
        input_m = build_fwd_binary_matrix2<F>(i, template_m);
      result_matrix_t fa = F::template apply<result_matrix_t>(input_1, 
                                                              input_m);
      EXPECT_EQ(input_m.size(), fa.size());
      expect_val_deriv_eq(F::apply_base(input_1, input_m(j)), fa(j));
    } 
  }   
}

template <typename F, typename FV, typename input_t1, typename input_t2, 
          int R, int C> 
void expect_fwd_binary_matrix_scalar_eq(
  const Eigen::Matrix<input_t1, R, C>& template_m, 
  const std::vector<input_t2>& template_scalar_v,
  bool seed_one = 1, bool seed_two = 1) {

  using std::vector;

  typedef Eigen::Matrix<input_t1, R, C> input_matrix_t;
  typedef Eigen::Matrix<FV, R, C> result_matrix_t;

  for (size_t i = 0; i < template_scalar_v.size(); ++i) {
    input_t2 input_2 = 0;
    if (seed_two)
      input_2 = build_binary_vector2<F>(template_scalar_v, i)[i];
    else
      input_2 = build_binary_vector2<F>(template_scalar_v)[i];
    for (int j = 0; j < template_m.size(); ++j) {
      input_matrix_t input_m; 
      if (seed_one)
        input_m = build_fwd_binary_matrix1<F>(i, template_m, j);
      else
        input_m = build_fwd_binary_matrix1<F>(i, template_m);
      result_matrix_t fa = F::template apply<result_matrix_t>(input_m, 
                                                              input_2);
      EXPECT_EQ(input_m.size(), fa.size());
      expect_val_deriv_eq(F::apply_base(input_m(j), input_2), fa(j));
    } 
  }   
}

template <typename F, typename FV, typename matrix_t1, typename matrix_t2>
void expect_fwd_binary_matrix_matrix_eq(const matrix_t1& template_m1, 
                                        const matrix_t2& template_m2, 
                                        bool seed_one = 1, 
                                        bool seed_two = 1) {

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
    result_matrix_t fa = F::template apply<result_matrix_t>(input_m1, 
                                                            input_m2);
    EXPECT_EQ(input_m1.size(), fa.size());
    EXPECT_EQ(input_m2.size(), fa.size());
    expect_val_deriv_eq(F::apply_base(input_m1(i), input_m2(i)), fa(i));
  }   

  if (template_m1.rows() > 1 && template_m1.cols() > 1) {
    matrix_t1 input_m1 = build_fwd_binary_matrix1<F>(template_m1);
    matrix_t2 input_m2 = build_fwd_binary_matrix2<F>(template_m2);
 
    result_matrix_t fb = F::template apply<result_matrix_t>(
      input_m1.block(1, 1, 1, 1), input_m2.block(1, 1, 1, 1)); 
    expect_val_deriv_eq(F::apply_base(input_m1(1, 1), input_m2(1, 1)), 
                        fb(0, 0));  
  }
}

template <typename F, typename FV, typename input_t1, typename input_t2, 
          int R, int C> void 
expect_fwd_binary_scalar_std_vector_matrix_eq(
  const std::vector<input_t1>& template_scalar_v, 
  const Eigen::Matrix<input_t2, R, C>& template_m, 
  bool seed_one = 1, bool seed_two = 1) {

  using std::vector;

  typedef Eigen::Matrix<input_t2, R, C> input_matrix_t;
  typedef Eigen::Matrix<FV, R, C> result_matrix_t;

  const size_t num_v = 2;
  for (size_t i = 0; i < template_scalar_v.size(); ++i) {
    input_t1 input_1 = 0; 
    if (seed_one)
      input_1 = build_binary_vector1<F>(template_scalar_v, i)[i];
    else
      input_1 = build_binary_vector1<F>(template_scalar_v)[i];
    for (size_t j = 0; j < num_v; ++j) {
      for (int k = 0; k < template_m.size(); ++k) {
        vector<input_matrix_t> input_mv; 
        for (size_t l = 0; l < num_v; ++l) {
          if (seed_two && l == j)
            input_mv.push_back(build_fwd_binary_matrix2<F>(
              i, template_m, k));
          else
            input_mv.push_back(build_fwd_binary_matrix2<F>(i, template_m));
        }
        vector<result_matrix_t> fa = F::template 
        apply<vector<result_matrix_t> >(input_1, input_mv);
        EXPECT_EQ(input_mv.size(), fa.size());
        EXPECT_EQ(input_mv[j].size(), fa[j].size());
        expect_val_deriv_eq(F::apply_base(input_1, input_mv[j](k)), 
                            fa[j](k));
      } 
    }   
  }
}

template <typename F, typename FV, typename input_t1, typename input_t2, 
          int R, int C> void 
expect_fwd_binary_std_vector_matrix_scalar_eq(
  const Eigen::Matrix<input_t1, R, C>& template_m, 
  const std::vector<input_t2>& template_scalar_v, 
  bool seed_one = 1, bool seed_two = 1) {

  using std::vector;

  typedef Eigen::Matrix<input_t1, R, C> input_matrix_t;
  typedef Eigen::Matrix<FV, R, C> result_matrix_t;

  const size_t num_v = 2;
  for (size_t i = 0; i < template_scalar_v.size(); ++i) {
    input_t2 input_2 = 0;
    if (seed_two)
      input_2 = build_binary_vector2<F>(template_scalar_v, i)[i];
    else
      input_2 = build_binary_vector2<F>(template_scalar_v)[i];
    for (size_t j = 0; j < num_v; ++j) {
      for (int k = 0; k < template_m.size(); ++k) {
        vector<input_matrix_t> input_mv;
        for (size_t l = 0; l < num_v; ++l) {
          if (seed_one && l == j)
            input_mv.push_back(build_fwd_binary_matrix1<F>(
              i, template_m, k));
          else
            input_mv.push_back(build_fwd_binary_matrix1<F>(i, template_m));
        }
        vector<result_matrix_t> fa = 
          F::template apply<vector<result_matrix_t> >(input_mv, input_2);
        EXPECT_EQ(input_mv.size(), fa.size());
        EXPECT_EQ(input_mv[j].size(), fa[j].size());
        expect_val_deriv_eq(F::apply_base(input_mv[j](k), input_2), 
                            fa[j](k));
      } 
    }   
  }
}

template <typename F, typename FV, typename matrix_t1, typename matrix_t2> 
void expect_fwd_binary_std_vector_matrix_std_vector_matrix_eq(
  const matrix_t1& template_m1, const matrix_t2& template_m2,
  bool seed_one = 1, bool seed_two = 1) {

  using std::vector;

  typedef Eigen::Matrix<FV, matrix_t1::RowsAtCompileTime, 
                        matrix_t1::ColsAtCompileTime> result_matrix_t;

  const size_t num_v = 2;
  for (size_t i = 0; i < num_v; ++i) {
    for (int j = 0; j < template_m1.size(); ++j) {
      vector<matrix_t1> input_mv1;
      for (size_t k = 0; k < num_v; ++k) {
        if (seed_one && k == i)
          input_mv1.push_back(build_fwd_binary_matrix1<F>(template_m1, j));
        else
          input_mv1.push_back(build_fwd_binary_matrix1<F>(template_m1));
      }
      vector<matrix_t2> input_mv2;
      for (size_t k = 0; k < num_v; ++k) {
        if (seed_one && k == i)
          input_mv2.push_back(build_fwd_binary_matrix2<F>(template_m2, j));
        else
          input_mv2.push_back(build_fwd_binary_matrix2<F>(template_m2));
      }
      vector<result_matrix_t> fa = F::template 
      apply<vector<result_matrix_t> >(input_mv1, input_mv2);
      EXPECT_EQ(input_mv1.size(), fa.size());
      EXPECT_EQ(input_mv2.size(), fa.size());
      EXPECT_EQ(input_mv1[i].size(), fa[i].size());
      EXPECT_EQ(input_mv2[i].size(), fa[i].size());
      expect_val_deriv_eq(F::apply_base(input_mv1[i](j), input_mv2[i](j)),
                          fa[i](j));
    }
  }   
}

template <typename F, typename FV, int R, int C> 
void expect_fwd_binary_scalar_matrix_all_eq(
  const std::vector<int>& int_template_v, 
  const std::vector<double>& d_template_v, 
  const std::vector<FV> fvar_template_v, 
  const Eigen::Matrix<double, R, C>& d_scalar_template_m, 
  const Eigen::Matrix<FV, R, C>& fvar_scalar_template_m) {

  expect_fwd_binary_scalar_matrix_eq<F, FV>(int_template_v, 
                                            fvar_scalar_template_m);  
  expect_fwd_binary_matrix_scalar_eq<F, FV>(fvar_scalar_template_m,
                                            int_template_v); 

  expect_fwd_binary_scalar_matrix_eq<F, FV>(fvar_template_v, 
                                            d_scalar_template_m);  
  expect_fwd_binary_matrix_scalar_eq<F, FV>(d_scalar_template_m, 
                                            fvar_template_v);  
  expect_fwd_binary_scalar_matrix_eq<F, FV>(d_template_v, 
                                            fvar_scalar_template_m);  
  expect_fwd_binary_matrix_scalar_eq<F, FV>(fvar_scalar_template_m, 
                                            d_template_v);

  expect_fwd_binary_scalar_matrix_eq<F, FV>(fvar_template_v, 
                                            fvar_scalar_template_m, 1, 0);  
  expect_fwd_binary_scalar_matrix_eq<F, FV>(fvar_template_v, 
                                            fvar_scalar_template_m, 0, 1);  
  expect_fwd_binary_scalar_matrix_eq<F, FV>(fvar_template_v, 
                                            fvar_scalar_template_m, 1, 1);  
  expect_fwd_binary_matrix_scalar_eq<F, FV>(fvar_scalar_template_m,
                                            fvar_template_v, 1, 0);  
  expect_fwd_binary_matrix_scalar_eq<F, FV>(fvar_scalar_template_m,
                                            fvar_template_v, 0, 1);  
  expect_fwd_binary_matrix_scalar_eq<F, FV>(fvar_scalar_template_m,
                                            fvar_template_v, 1, 1);  
}

template <typename F, typename FV, int R, int C> void 
expect_fwd_binary_scalar_std_vector_matrix_all_eq(
  const std::vector<int>& int_template_v, 
  const std::vector<double>& d_template_v,
  const std::vector<FV> fvar_template_v, 
  const Eigen::Matrix<double, R, C>& d_scalar_template_m, 
  const Eigen::Matrix<FV, R, C>& fvar_scalar_template_m) {

  expect_fwd_binary_scalar_std_vector_matrix_eq<F, FV>(
    int_template_v, fvar_scalar_template_m);
  expect_fwd_binary_std_vector_matrix_scalar_eq<F, FV>(
    fvar_scalar_template_m, int_template_v);

  expect_fwd_binary_scalar_std_vector_matrix_eq<F, FV>(
    fvar_template_v, d_scalar_template_m);
  expect_fwd_binary_std_vector_matrix_scalar_eq<F, FV>(
    d_scalar_template_m, fvar_template_v);
  expect_fwd_binary_scalar_std_vector_matrix_eq<F, FV>(
    d_template_v, fvar_scalar_template_m);
  expect_fwd_binary_std_vector_matrix_scalar_eq<F, FV>(
    fvar_scalar_template_m, d_template_v);

  expect_fwd_binary_scalar_std_vector_matrix_eq<F, FV>(
    fvar_template_v, fvar_scalar_template_m, 1, 0);
  expect_fwd_binary_scalar_std_vector_matrix_eq<F, FV>(
    fvar_template_v, fvar_scalar_template_m, 0, 1);
  expect_fwd_binary_scalar_std_vector_matrix_eq<F, FV>(
    fvar_template_v, fvar_scalar_template_m, 1, 1);
  expect_fwd_binary_std_vector_matrix_scalar_eq<F, FV>(
    fvar_scalar_template_m, fvar_template_v, 1, 0);
  expect_fwd_binary_std_vector_matrix_scalar_eq<F, FV>(
    fvar_scalar_template_m, fvar_template_v, 0, 1);
  expect_fwd_binary_std_vector_matrix_scalar_eq<F, FV>(
    fvar_scalar_template_m, fvar_template_v, 1, 1);
}

template <typename F, typename FV, int R, int C> void 
expect_fwd_binary_matrix_value(Eigen::Matrix<double, R, C> template_m) {
  using std::vector;
  
  vector<int> int_template_v;
  vector<double> d_template_v;
  vector<FV> fvar_template_v;

  typedef Eigen::Matrix<double, R, C> d_matrix_t;
  typedef Eigen::Matrix<FV, R, C> fvar_matrix_t;

  d_matrix_t d_scalar_template_m;
  d_matrix_t d_template_m;
  fvar_matrix_t fvar_scalar_template_m;
  fvar_matrix_t fvar_template_m;

  d_scalar_template_m = build_template_matrix(d_scalar_template_m, 3, 5);
  d_template_m = build_template_matrix(d_scalar_template_m,
                                       3, F::valid_inputs1().size());
  fvar_scalar_template_m = build_template_matrix(fvar_scalar_template_m,
                                                 3, 5);
  fvar_template_m = build_template_matrix(fvar_scalar_template_m,
                                          3, F::valid_inputs1().size());

  expect_fwd_binary_scalar_matrix_all_eq<F, FV>(int_template_v, 
                                                d_template_v, 
                                                fvar_template_v, 
                                                d_scalar_template_m, 
                                                fvar_scalar_template_m);  

  expect_fwd_binary_matrix_matrix_eq<F, FV>(fvar_template_m, d_template_m);
  expect_fwd_binary_matrix_matrix_eq<F, FV>(d_template_m, fvar_template_m);
  expect_fwd_binary_matrix_matrix_eq<F, FV>(fvar_template_m, 
                                            fvar_template_m, 1, 0);
  expect_fwd_binary_matrix_matrix_eq<F, FV>(fvar_template_m,
                                            fvar_template_m, 0, 1);
  expect_fwd_binary_matrix_matrix_eq<F, FV>(fvar_template_m, 
                                            fvar_template_m, 1, 1);

  expect_fwd_binary_scalar_std_vector_matrix_all_eq<F, FV>(
    int_template_v, d_template_v, fvar_template_v, d_scalar_template_m, 
    fvar_scalar_template_m);

  expect_fwd_binary_std_vector_matrix_std_vector_matrix_eq<F, FV>(
    fvar_template_m, d_template_m);
  expect_fwd_binary_std_vector_matrix_std_vector_matrix_eq<F, FV>(
    d_template_m, fvar_template_m);
  expect_fwd_binary_std_vector_matrix_std_vector_matrix_eq<F, FV>(
    fvar_template_m, fvar_template_m, 1, 0);
  expect_fwd_binary_std_vector_matrix_std_vector_matrix_eq<F, FV>(
    fvar_template_m, fvar_template_m, 0, 1);
  expect_fwd_binary_std_vector_matrix_std_vector_matrix_eq<F, FV>(
    fvar_template_m, fvar_template_m, 1, 1);
}

#endif
