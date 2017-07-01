#ifndef TEST_UNIT_MATH_FWD_MAT_VECTORIZE_EXPECT_FWD_BINARY_MATRIX_VALUE_HPP
#define TEST_UNIT_MATH_FWD_MAT_VECTORIZE_EXPECT_FWD_BINARY_MATRIX_VALUE_HPP

#include <stan/math/fwd/core/fvar.hpp>
#include <test/unit/math/fwd/mat/vectorize/expect_fwd_binary_scalar_matrix_all_eq.hpp>
#include <test/unit/math/fwd/mat/vectorize/expect_fwd_binary_matrix_matrix_eq.hpp>
#include <test/unit/math/fwd/mat/vectorize/expect_fwd_binary_scalar_std_vector_matrix_all_eq.hpp>
#include <test/unit/math/fwd/mat/vectorize/expect_fwd_binary_std_vector_matrix_std_vector_matrix_eq.hpp>
#include <Eigen/Dense>
#include <vector>

template <typename F, typename FV>
void expect_fwd_binary_matrix_value() {
  using std::vector;
  
  vector<int> int_template_v;
  vector<double> d_template_v;
  vector<FV> fvar_template_v;

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> 
  d_scalar_template_m(3, 5);
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> 
  d_template_m(3, F::valid_inputs1().size());
  Eigen::Matrix<FV, Eigen::Dynamic, Eigen::Dynamic> 
  fvar_scalar_template_m(3, 5);
  Eigen::Matrix<FV, Eigen::Dynamic, Eigen::Dynamic> 
  fvar_template_m(3, F::valid_inputs1().size());

  expect_fwd_binary_scalar_matrix_all_eq<F, FV>(int_template_v, 
  d_template_v, fvar_template_v, d_scalar_template_m, 
  fvar_scalar_template_m);  

  expect_fwd_binary_matrix_matrix_eq<F, FV>(fvar_template_m, d_template_m);
  expect_fwd_binary_matrix_matrix_eq<F, FV>(d_template_m, fvar_template_m);
  expect_fwd_binary_matrix_matrix_eq<F, FV>(fvar_template_m, 
  fvar_template_m, 1, 0);
  expect_fwd_binary_matrix_matrix_eq<F, FV>(fvar_template_m,
  fvar_template_m, 0, 1);
  expect_fwd_binary_matrix_matrix_eq<F, FV>(fvar_template_m, 
  fvar_template_m, 1, 1);

  expect_fwd_binary_scalar_std_vector_matrix_all_eq<F, FV>(int_template_v,
  d_template_v, fvar_template_v, d_scalar_template_m, 
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
