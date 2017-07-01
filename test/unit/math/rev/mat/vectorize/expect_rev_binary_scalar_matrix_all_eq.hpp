#ifndef TEST_UNIT_MATH_REV_MAT_VECTORIZE_EXPECT_REV_BINARY_SCALAR_MATRIX_ALL_EQ
#define TEST_UNIT_MATH_REV_MAT_VECTORIZE_EXPECT_REV_BINARY_SCALAR_MATRIX_ALL_EQ

#include <test/unit/math/rev/mat/vectorize/expect_rev_binary_scalar_matrix_eq.hpp>
#include <test/unit/math/rev/mat/vectorize/expect_rev_binary_matrix_scalar_eq.hpp>
#include <Eigen/Dense>
#include <vector>

template <typename F, int R, int C> 
void expect_rev_binary_scalar_matrix_all_eq(
const std::vector<int>& int_template_v, const std::vector<double>&
d_template_v, const std::vector<stan::math::var>& var_template_v,
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
#endif
