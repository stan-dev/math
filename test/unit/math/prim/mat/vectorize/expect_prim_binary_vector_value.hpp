#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_PRIM_BINARY_VECTOR_VALUE_HPP
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_PRIM_BINARY_VECTOR_VALUE_HPP

#include <test/unit/math/prim/mat/vectorize/expect_prim_binary_scalar_matrix_eq.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_prim_binary_matrix_scalar_eq.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_prim_binary_matrix_matrix_eq.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_prim_binary_scalar_std_vector_matrix_eq.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_prim_binary_std_vector_matrix_scalar_eq.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_prim_binary_std_vector_matrix_std_vector_matrix_eq.hpp>
#include <vector>
#include <Eigen/Dense>

template <typename F, typename vector_t>
void expect_prim_binary_vector_value() {
  using std::vector;
  using Eigen::MatrixXd;

  std::vector<int> int_valid_inputs1 = F::int_valid_inputs1();
  std::vector<int> int_valid_inputs2 = F::int_valid_inputs2();
  std::vector<double> valid_inputs1 = F::valid_inputs1();
  std::vector<double> valid_inputs2 = F::valid_inputs2();

  //Tests int, vector
  for (size_t i = 0; i < int_valid_inputs1.size(); ++i) {
    vector_t b1 = vector_t::Constant(5, valid_inputs1[i]);
    vector_t b2 = vector_t::Constant(5, valid_inputs2[i]);
    expect_prim_binary_scalar_matrix_eq<F>(int_valid_inputs1[i], b2);
    expect_prim_binary_matrix_scalar_eq<F>(b1, int_valid_inputs2[i]);
  }
  
  //Tests double, vector
  for (size_t i = 0; i < valid_inputs1.size(); ++i) {
    vector_t b1 = vector_t::Constant(5, valid_inputs1[i]);
    vector_t b2 = vector_t::Constant(5, valid_inputs2[i]);
    expect_prim_binary_scalar_matrix_eq<F>(valid_inputs1[i], b2);
    expect_prim_binary_matrix_scalar_eq<F>(b1, valid_inputs2[i]);
  }
  
  //Tests vector, vector
  vector_t b1 = vector_t::Map(valid_inputs1.data(), valid_inputs1.size());
  vector_t b2 = vector_t::Map(valid_inputs2.data(), valid_inputs2.size());
  expect_prim_binary_matrix_matrix_eq<F>(b1, b2);

  //Tests int, std::vector<vector>
  for (size_t i = 0; i < int_valid_inputs1.size(); ++i) {
    vector_t c1 = vector_t::Map(valid_inputs1.data(), valid_inputs1.size());
    vector_t c2 = vector_t::Map(valid_inputs2.data(), valid_inputs2.size());
    vector<vector_t> d1;
    d1.push_back(c1);
    d1.push_back(c1);
    vector<vector_t> d2;
    d2.push_back(c2);
    d2.push_back(c2);
    expect_prim_binary_scalar_std_vector_matrix_eq<F>(
      int_valid_inputs1[i], d2);
    expect_prim_binary_std_vector_matrix_scalar_eq<F>(
      d1, int_valid_inputs2[i]);
  }

  //Tests double, std::vector<vector>
  for (size_t i = 0; i < valid_inputs1.size(); ++i) {
    vector_t c1 = vector_t::Map(valid_inputs1.data(), valid_inputs1.size());
    vector_t c2 = vector_t::Map(valid_inputs2.data(), valid_inputs2.size());
    vector<vector_t> d1;
    d1.push_back(c1);
    d1.push_back(c1);
    vector<vector_t> d2;
    d2.push_back(c2);
    d2.push_back(c2);
    expect_prim_binary_scalar_std_vector_matrix_eq<F>(
      valid_inputs1[i], d2);
    expect_prim_binary_std_vector_matrix_scalar_eq<F>(
      d1, valid_inputs2[i]);
  }

  //Tests std::vector<vector>, std::vector<vector>
  vector<vector_t> d1;
  d1.push_back(b1);
  d1.push_back(b1);
  vector<vector_t> d2;
  d2.push_back(b2);
  d2.push_back(b2);
  expect_prim_binary_std_vector_matrix_std_vector_matrix_eq<F>(d1, d2);
}
#endif
