#ifndef TEST_UNIT_MATH_FWD_MAT_VECTORIZE_EXPECT_FWD_BINARY_MATRIX_VALUE_HPP
#define TEST_UNIT_MATH_FWD_MAT_VECTORIZE_EXPECT_FWD_BINARY_MATRIX_VALUE_HPP

#include <stan/math/fwd/core/fvar.hpp>
#include <test/unit/math/fwd/mat/vectorize/expect_fwd_binary_scalar_matrix_eq.hpp>
#include <test/unit/math/fwd/mat/vectorize/expect_fwd_binary_matrix_scalar_eq.hpp>
#include <test/unit/math/fwd/mat/vectorize/expect_fwd_binary_matrix_matrix_eq.hpp>
#include <test/unit/math/fwd/mat/vectorize/expect_fwd_binary_scalar_std_vector_matrix_eq.hpp>
#include <test/unit/math/fwd/mat/vectorize/expect_fwd_binary_std_vector_matrix_scalar_eq.hpp>
#include <test/unit/math/fwd/mat/vectorize/expect_fwd_binary_std_vector_matrix_std_vector_matrix_eq.hpp>
#include <Eigen/Dense>
#include <vector>

template <typename F, typename FV>
void expect_fwd_binary_matrix_value() {
  using std::vector;
  
  vector<int> int_template_v(F::int_valid_inputs1().size());
  vector<double> d_template_v(F::valid_inputs1().size());
  vector<FV> fvar_template_v(F::valid_inputs1().size());

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> 
  d_scalar_template_m(3, 5);
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> 
  d_template_m(3, F::valid_inputs1().size());
  Eigen::Matrix<FV, Eigen::Dynamic, Eigen::Dynamic> 
  fvar_scalar_template_m(3, 5);
  Eigen::Matrix<FV, Eigen::Dynamic, Eigen::Dynamic> 
  fvar_template_m(3, F::valid_inputs1().size());

  //scalar, matrix
  //fvar, int
  expect_fwd_binary_scalar_matrix_eq<F, FV>(int_template_v, 
  fvar_scalar_template_m);  
  expect_fwd_binary_matrix_scalar_eq<F, FV>(fvar_scalar_template_m,
  int_template_v); 
  //fvar, double 
  expect_fwd_binary_scalar_matrix_eq<F, FV>(fvar_template_v, 
  d_scalar_template_m);  
  expect_fwd_binary_matrix_scalar_eq<F, FV>(d_scalar_template_m, 
  fvar_template_v);  
  expect_fwd_binary_scalar_matrix_eq<F, FV>(d_template_v, 
  fvar_scalar_template_m);  
  expect_fwd_binary_matrix_scalar_eq<F, FV>(fvar_scalar_template_m, 
  d_template_v);
  //fvar, fvar  
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

  //matrix, matrix
  expect_fwd_binary_matrix_matrix_eq<F, FV>(fvar_template_m, d_template_m);
  expect_fwd_binary_matrix_matrix_eq<F, FV>(d_template_m, fvar_template_m);
  expect_fwd_binary_matrix_matrix_eq<F, FV>(fvar_template_m, 
  fvar_template_m, 1, 0);
  expect_fwd_binary_matrix_matrix_eq<F, FV>(fvar_template_m,
  fvar_template_m, 0, 1);
  expect_fwd_binary_matrix_matrix_eq<F, FV>(fvar_template_m, 
  fvar_template_m, 1, 1);

  //scalar, matrix
  //fvar, int
  expect_fwd_binary_scalar_std_vector_matrix_eq<F, FV>(int_template_v,
  fvar_scalar_template_m);
  expect_fwd_binary_std_vector_matrix_scalar_eq<F, FV>(
  fvar_scalar_template_m, int_template_v);
  //fvar, double
  expect_fwd_binary_scalar_std_vector_matrix_eq<F, FV>(fvar_template_v,
  d_scalar_template_m);
  expect_fwd_binary_std_vector_matrix_scalar_eq<F, FV>(d_scalar_template_m, 
  fvar_template_v);
  expect_fwd_binary_scalar_std_vector_matrix_eq<F, FV>(d_template_v,
  fvar_scalar_template_m);
  expect_fwd_binary_std_vector_matrix_scalar_eq<F, FV>(
  fvar_scalar_template_m, d_template_v);
  //fvar, fvar
  expect_fwd_binary_scalar_std_vector_matrix_eq<F, FV>(fvar_template_v,
  fvar_scalar_template_m, 1, 0);
  expect_fwd_binary_scalar_std_vector_matrix_eq<F, FV>(fvar_template_v,
  fvar_scalar_template_m, 0, 1);
  expect_fwd_binary_scalar_std_vector_matrix_eq<F, FV>(fvar_template_v,
  fvar_scalar_template_m, 1, 1);
  expect_fwd_binary_std_vector_matrix_scalar_eq<F, FV>(
  fvar_scalar_template_m, fvar_template_v, 1, 0);
  expect_fwd_binary_std_vector_matrix_scalar_eq<F, FV>(
  fvar_scalar_template_m, fvar_template_v, 0, 1);
  expect_fwd_binary_std_vector_matrix_scalar_eq<F, FV>(
  fvar_scalar_template_m, fvar_template_v, 1, 1);

  //vector<matrix>, vector<matrix>
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
