#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_PRIM_BINARY_STD_VECTOR_VALUE_HPP
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_PRIM_BINARY_STD_VECTOR_VALUE_HPP

#include <test/unit/math/prim/mat/vectorize/expect_val_eq.hpp>
#include <vector>
#include <gtest/gtest.h>

template <typename F, typename vector_t1, typename vector_t2>
void expect_binary_std_vector_std_vector_eq(
    std::vector<vector_t1> input_v1, std::vector<vector_t2> input_v2) {
  using std::vector;

  vector<double> fc = F::template apply<vector<double> >(input_v1,
                                                         input_v2);
  EXPECT_EQ(input_v1.size(), fc.size());
  EXPECT_EQ(input_v2.size(), fc.size());
  for (size_t i = 0; i < input_v1.size(); ++i) {
    double exp_v = F::apply_base(input_v1[i], input_v2[i]);
    expect_val_eq(exp_v, fc[i]);
  }
}

template <typename F, typename vector_t1, typename vector_t2>
void expect_binary_std_vector_std_vector_std_vector_std_vector_eq(
    std::vector<std::vector<vector_t1> > input_v1,
    std::vector<std::vector<vector_t2> > input_v2) {
  using std::vector;

  vector<vector<double> > fd = F::template apply<vector<vector<double> > >(
                                   input_v1, input_v2);
  EXPECT_EQ(input_v1.size(), fd.size());
  EXPECT_EQ(input_v2.size(), fd.size());
  for (size_t i = 0; i < input_v1.size(); ++i) {
    EXPECT_EQ(input_v1[i].size(), fd[i].size());
    EXPECT_EQ(input_v2[i].size(), fd[i].size());
    for (size_t j = 0; j < input_v1[i].size(); ++j) {
      double exp_v = F::apply_base(input_v1[i][j], input_v2[i][j]);
      expect_val_eq(exp_v, fd[i][j]);
    }
  }
}

template <typename F>
void expect_prim_binary_std_vector_value() {
  using std::vector;

  vector<int> int_valid_inputs1 = F::int_valid_inputs1();
  vector<int> int_valid_inputs2 = F::int_valid_inputs2();
  vector<double> valid_inputs1 = F::valid_inputs1();
  vector<double> valid_inputs2 = F::valid_inputs2();

  expect_binary_std_vector_std_vector_eq<F>(int_valid_inputs1,
                                            valid_inputs2);
  expect_binary_std_vector_std_vector_eq<F>(valid_inputs1,
                                            int_valid_inputs2);
  expect_binary_std_vector_std_vector_eq<F>(int_valid_inputs1,
                                            int_valid_inputs2);
  expect_binary_std_vector_std_vector_eq<F>(valid_inputs1, valid_inputs2);

  vector<vector<int> > a1;
  a1.push_back(int_valid_inputs1);
  a1.push_back(int_valid_inputs1);
  vector<vector<int> > a2;
  a2.push_back(int_valid_inputs2);
  a2.push_back(int_valid_inputs2);
  vector<vector<double> > b1;
  b1.push_back(valid_inputs1);
  b1.push_back(valid_inputs1);
  vector<vector<double> > b2;
  b2.push_back(valid_inputs2);
  b2.push_back(valid_inputs2);

  //vector<vector>, vector<vector> tests
  expect_binary_std_vector_std_vector_std_vector_std_vector_eq<F>(a1, b2);
  expect_binary_std_vector_std_vector_std_vector_std_vector_eq<F>(b1, a2);
  expect_binary_std_vector_std_vector_std_vector_std_vector_eq<F>(a1, a2);
  expect_binary_std_vector_std_vector_std_vector_std_vector_eq<F>(b2, b2);
}

#endif
