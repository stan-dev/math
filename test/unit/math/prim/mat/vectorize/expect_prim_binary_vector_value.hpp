#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_PRIM_BINARY_VECTOR_VALUE_HPP
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_PRIM_BINARY_VECTOR_VALUE_HPP

#include <test/unit/math/prim/mat/vectorize/expect_prim_binary_scalar_matrix_eq.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_prim_binary_matrix_matrix_eq.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_prim_binary_scalar_std_vector_matrix_eq.hpp>
#include <test/unit/math/prim/mat/vectorize/expect_prim_binary_std_vector_matrix_std_vector_matrix_eq.hpp>
#include <vector>
#include <Eigen/Dense>

template <typename F, typename vector_t>
void expect_prim_binary_vector_value() {
  using std::vector;

  std::vector<int> int_valid_inputs = F::int_valid_inputs();
  std::vector<double> valid_inputs = F::valid_inputs();

  vector_t b = vector_t::Map(valid_inputs.data(), valid_inputs.size());
  
  //Tests int, vector
  for (size_t i = 0; i < int_valid_inputs.size(); ++i) {
    expect_prim_binary_scalar_matrix_eq<F>(int_valid_inputs[i], b);
  }
  
  //Tests double, vector
  for (size_t i = 0; i < valid_inputs.size(); ++i) {
    expect_prim_binary_scalar_matrix_eq<F>(valid_inputs[i], b);
  }
  
  //Tests vector, vector
  expect_prim_binary_matrix_matrix_eq<F>(b);

  //Tests int, std::vector<vector>
  vector<vector_t> c;
  c.push_back(b);
  c.push_back(b);
  for (size_t i = 0; i < int_valid_inputs.size(); ++i) {
    expect_prim_binary_scalar_std_vector_matrix_eq<F>(
    int_valid_inputs[i], c, b.size());
  }

  //Tests double, std::vector<vector>
  for (size_t i = 0; i < valid_inputs.size(); ++i) {
    expect_prim_binary_scalar_std_vector_matrix_eq<F>(
    valid_inputs[i], c, b.size());
  }

  //Tests std::vector<vector>, std::vector<vector>
  expect_prim_binary_std_vector_matrix_std_vector_matrix_eq<F>(c, b.size());
}
#endif
