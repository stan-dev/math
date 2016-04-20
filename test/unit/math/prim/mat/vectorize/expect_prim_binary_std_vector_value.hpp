#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_PRIM_BINARY_STD_VECTOR_VALUE_HPP
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_PRIM_BINARY_STD_VECTOR_VALUE_HPP

#include <stan/math/prim/scal/fun/is_nan.hpp>
#include <vector>
#include <gtest/gtest.h>

template <typename F, typename input_t, typename vector_t> 
void expect_binary_scalar_std_vector_eq(input_t input, 
                                        std::vector<vector_t> input_v) {
  using stan::math::is_nan;
  using std::vector;

  vector<double> fa
  = F::template apply<vector<double> >(input, input_v); 
  EXPECT_EQ(input_v.size(), fa.size()); 
  for (size_t i = 0; i < fa.size(); ++i) {
    double exp_v = F::apply_base(input, input_v[i]);
    if (is_nan(exp_v) && is_nan(fa[i])) continue;
    EXPECT_FLOAT_EQ(exp_v, fa[i]);  
  }

  fa = F::template apply<vector<double> >(input_v, input);
  EXPECT_EQ(input_v.size(), fa.size()); 
  for (size_t i = 0; i < fa.size(); ++i) {
    double exp_v = F::apply_base(input_v[i], input);
    if (is_nan(exp_v) && is_nan(fa[i])) continue;
    EXPECT_FLOAT_EQ(exp_v, fa[i]);  
  }
}

template <typename F, typename input_t, typename vector_t> 
void expect_binary_scalar_nested_std_vector_eq(
input_t input, std::vector<std::vector<vector_t> > input_v, size_t v_size) {
  using stan::math::is_nan;
  using std::vector;

  vector<vector<double> > fb
  = F::template apply<vector<vector<double> > >(input, input_v); 
  EXPECT_EQ(input_v.size(), fb.size()); 
  for (size_t i = 0; i < fb.size(); ++i) {
    EXPECT_EQ(v_size, fb[i].size());
    for (size_t j = 0; j < fb[i].size(); ++j) {
      double exp_v = F::apply_base(input, input_v[i][j]);
      if (is_nan(exp_v) && is_nan(fb[i][j])) continue;
      EXPECT_FLOAT_EQ(exp_v, fb[i][j]);  
    }
  }

  fb = F::template apply<vector<vector<double> > >(input_v, input); 
  EXPECT_EQ(input_v.size(), fb.size()); 
  for (size_t i = 0; i < fb.size(); ++i) {
    EXPECT_EQ(v_size, fb[i].size());
    for (size_t j = 0; j < fb[i].size(); ++j) {
      double exp_v = F::apply_base(input_v[i][j], input);
      if (is_nan(exp_v) && is_nan(fb[i][j])) continue;
      EXPECT_FLOAT_EQ(exp_v, fb[i][j]);  
    }
  }
}

template <typename F, typename vector_t1, typename vector_t2> 
void expect_binary_std_vector_std_vector_eq(
std::vector<vector_t1> input_v1, std::vector<vector_t2> input_v2) {
  using stan::math::is_nan;
  using std::vector;

  vector<double> fc = F::template apply<vector<double> >
                      (input_v1, input_v2);
  EXPECT_EQ(input_v1.size(), fc.size());
  EXPECT_EQ(input_v2.size(), fc.size());
  for (size_t i = 0; i < fc.size(); ++i) {
    double exp_v = F::apply_base(input_v1[i], input_v2[i]);
    if (is_nan(exp_v) && is_nan(fc[i])) continue;
    EXPECT_FLOAT_EQ(exp_v, fc[i]);
  }
}

template <typename F, typename vector_t1, typename vector_t2> 
void expect_binary_nested_std_vector_nested_std_vector_eq(
std::vector<std::vector<vector_t1> > input_v1, std::vector<
std::vector<vector_t2> > input_v2, size_t v1_size, size_t v2_size) {
  using stan::math::is_nan;
  using std::vector;

  vector<vector<double> > fd = F::template apply<vector<vector<double> > >
                               (input_v1, input_v2);
  EXPECT_EQ(input_v1.size(), fd.size());
  EXPECT_EQ(input_v2.size(), fd.size());
  for (size_t i = 0; i < fd.size(); ++i) {
    EXPECT_EQ(v1_size, fd[i].size());
    EXPECT_EQ(v2_size, fd[i].size());
    for (size_t j = 0; j < fd[i].size(); ++j) {
      double exp_v = F::apply_base(input_v1[i][j], input_v2[i][j]);
      if (is_nan(exp_v) && is_nan(fd[i][j])) continue;
      EXPECT_FLOAT_EQ(exp_v, fd[i][j]);
    }
  }
}

template <typename F>
void expect_prim_binary_std_vector_value() {
  using stan::math::is_nan;
  using std::vector;

  vector<int> int_valid_inputs = F::int_valid_inputs();
  vector<double> valid_inputs = F::valid_inputs();
  //int, vector tests
  for (size_t i = 0; i < valid_inputs.size(); ++i) {
    int input = int_valid_inputs[i];
    expect_binary_scalar_std_vector_eq<F>(input, int_valid_inputs);
    expect_binary_scalar_std_vector_eq<F>(input, valid_inputs);
  }

  //double, vector tests
  for (size_t i = 0; i < valid_inputs.size(); ++i) {
    double input = valid_inputs[i];
    expect_binary_scalar_std_vector_eq<F>(input, int_valid_inputs);
    expect_binary_scalar_std_vector_eq<F>(input, valid_inputs);
  }
  //vector, vector tests
  expect_binary_std_vector_std_vector_eq<F>(int_valid_inputs, valid_inputs);
  expect_binary_std_vector_std_vector_eq<F>(valid_inputs, int_valid_inputs);
  expect_binary_std_vector_std_vector_eq<F>(int_valid_inputs, 
                                            int_valid_inputs);
  expect_binary_std_vector_std_vector_eq<F>(valid_inputs, valid_inputs);


  vector<vector<int> > a;
  a.push_back(int_valid_inputs);
  a.push_back(int_valid_inputs);
  vector<vector<double> > b;
  b.push_back(valid_inputs);
  b.push_back(valid_inputs);

  //int, vector<vector> tests
  for (size_t i = 0; i < int_valid_inputs.size(); ++i) {
    int input = int_valid_inputs[i];
    expect_binary_scalar_nested_std_vector_eq<F>(input, a,
                                                 int_valid_inputs.size());
    expect_binary_scalar_nested_std_vector_eq<F>(input, b, 
                                                 valid_inputs.size());
  }

  //double, vector<vector> tests
  for (size_t i = 0; i < valid_inputs.size(); ++i) {
    double input = valid_inputs[i];
    expect_binary_scalar_nested_std_vector_eq<F>(input, a,
                                                 int_valid_inputs.size());
    expect_binary_scalar_nested_std_vector_eq<F>(input, b, 
                                                 valid_inputs.size());
  }

  //vector<vector>, vector<vector> tests
  expect_binary_nested_std_vector_nested_std_vector_eq<F>(
  a, b, int_valid_inputs.size(), valid_inputs.size());
  expect_binary_nested_std_vector_nested_std_vector_eq<F>(
  b, a, valid_inputs.size(), int_valid_inputs.size());
  expect_binary_nested_std_vector_nested_std_vector_eq<F>(
  a, a, int_valid_inputs.size(), int_valid_inputs.size());
  expect_binary_nested_std_vector_nested_std_vector_eq<F>(
  b, b, valid_inputs.size(), valid_inputs.size());
}

#endif
