#ifndef TEST_UNIT_MATH_REV_MAT_VECTORIZE_EXPECT_REV_BINARY_MATRIX_VALUE_HPP
#define TEST_UNIT_MATH_REV_MAT_VECTORIZE_EXPECT_REV_BINARY_MATRIX_VALUE_HPP

#include <stan/math/rev/core/var.hpp>
#include <test/unit/math/rev/mat/vectorize/expect_rev_binary_scalar_matrix_eq.hpp>
#include <test/unit/math/rev/mat/vectorize/expect_rev_binary_matrix_scalar_eq.hpp>
#include <test/unit/math/rev/mat/vectorize/expect_rev_binary_matrix_matrix_eq.hpp>
#include <test/unit/math/rev/mat/vectorize/expect_rev_binary_scalar_std_vector_matrix_eq.hpp>
#include <test/unit/math/rev/mat/vectorize/expect_rev_binary_std_vector_matrix_scalar_eq.hpp>
#include <test/unit/math/rev/mat/vectorize/expect_rev_binary_std_vector_matrix_std_vector_matrix_eq.hpp>
#include <Eigen/Dense>
#include <vector>

template <typename F>
void expect_rev_binary_matrix_value() {
  using stan::math::var;
  using std::vector;
  
  vector<int> int_template_v(F::int_valid_inputs1().size());
  vector<double> d_template_v(F::valid_inputs1().size());
  vector<var> var_template_v(F::valid_inputs1().size());

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> 
  d_scalar_template_m(3, 5);
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> 
  d_template_m(3, F::valid_inputs1().size());
  Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> 
  var_scalar_template_m(3, 5);
  Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> 
  var_template_m(3, F::valid_inputs1().size());

  //scalar, matrix
  //var, int
  expect_rev_binary_scalar_matrix_eq<F>(int_template_v, 
  var_scalar_template_m);  
  expect_rev_binary_matrix_scalar_eq<F>(var_scalar_template_m,
  int_template_v); 
  //var, double 
  expect_rev_binary_scalar_matrix_eq<F>(var_template_v, 
  d_scalar_template_m);  
  expect_rev_binary_matrix_scalar_eq<F>(d_scalar_template_m, 
  var_template_v);  
  expect_rev_binary_scalar_matrix_eq<F>(d_template_v, 
  var_scalar_template_m);  
  expect_rev_binary_matrix_scalar_eq<F>(var_scalar_template_m, 
  d_template_v);
  //var, var  
  expect_rev_binary_scalar_matrix_eq<F>(var_template_v, 
  var_scalar_template_m);  
  expect_rev_binary_matrix_scalar_eq<F>(var_scalar_template_m,
  var_template_v);  

  //matrix, matrix
  expect_rev_binary_matrix_matrix_eq<F>(var_template_m, d_template_m);
  expect_rev_binary_matrix_matrix_eq<F>(d_template_m, var_template_m);
  expect_rev_binary_matrix_matrix_eq<F>(var_template_m, var_template_m);

  //scalar, matrix
  //var, int
  expect_rev_binary_scalar_std_vector_matrix_eq<F>(int_template_v,
  var_scalar_template_m);
  expect_rev_binary_std_vector_matrix_scalar_eq<F>(var_scalar_template_m, 
  int_template_v);
  //var, double
  expect_rev_binary_scalar_std_vector_matrix_eq<F>(var_template_v,
  d_scalar_template_m);
  expect_rev_binary_std_vector_matrix_scalar_eq<F>(d_scalar_template_m, 
  var_template_v);
  expect_rev_binary_scalar_std_vector_matrix_eq<F>(d_template_v,
  var_scalar_template_m);
  expect_rev_binary_std_vector_matrix_scalar_eq<F>(var_scalar_template_m, 
  d_template_v);
  //var, var
  expect_rev_binary_scalar_std_vector_matrix_eq<F>(var_template_v,
  var_scalar_template_m);
  expect_rev_binary_std_vector_matrix_scalar_eq<F>(var_scalar_template_m, 
  var_template_v);

  //vector<matrix>, vector<matrix>
  expect_rev_binary_std_vector_matrix_std_vector_matrix_eq<F>(
  var_template_m, d_template_m);
  expect_rev_binary_std_vector_matrix_std_vector_matrix_eq<F>(
  d_template_m, var_template_m);
  expect_rev_binary_std_vector_matrix_std_vector_matrix_eq<F>(
  var_template_m, var_template_m);
}

#endif
