#ifndef TEST_UNIT_MATH_FWD_MAT_VECTORIZE_EXPECT_FWD_BINARY_VECTOR_VALUE_HPP
#define TEST_UNIT_MATH_FWD_MAT_VECTORIZE_EXPECT_FWD_BINARY_VECTOR_VALUE_HPP

#include <stan/math/fwd/core/fvar.hpp>
#include <test/unit/math/fwd/mat/vectorize/expect_fwd_binary_scalar_matrix_eq.hpp>
#include <test/unit/math/fwd/mat/vectorize/expect_fwd_binary_matrix_scalar_eq.hpp>
#include <test/unit/math/fwd/mat/vectorize/expect_fwd_binary_matrix_matrix_eq.hpp>
#include <test/unit/math/fwd/mat/vectorize/expect_fwd_binary_scalar_std_vector_matrix_eq.hpp>
#include <test/unit/math/fwd/mat/vectorize/expect_fwd_binary_std_vector_matrix_scalar_eq.hpp>
#include <test/unit/math/fwd/mat/vectorize/expect_fwd_binary_std_vector_matrix_std_vector_matrix_eq.hpp>
#include <Eigen/Dense>
#include <vector>

template <typename F, typename FV, typename vector_t>
void expect_fwd_binary_vector_value() {
  using stan::math::fvar;
  using std::vector;
  
  vector<int> int_template_sv(F::int_valid_inputs1().size());
  vector<double> d_template_sv(F::valid_inputs1().size());
  vector<FV> fvar_template_sv(F::valid_inputs1().size());

  vector_t d_scalar_template_mv(5);
  vector_t d_template_mv(F::valid_inputs1().size());
  Eigen::Matrix<FV, vector_t::RowsAtCompileTime, 
  vector_t::ColsAtCompileTime> fvar_scalar_template_mv(5);
  Eigen::Matrix<FV, vector_t::RowsAtCompileTime, 
  vector_t::ColsAtCompileTime> fvar_template_mv(F::valid_inputs1().size());

  //scalar, vector
  //fvar, int
  expect_fwd_binary_scalar_matrix_eq<F, FV>(int_template_sv, 
  fvar_scalar_template_mv);  
  expect_fwd_binary_matrix_scalar_eq<F, FV>(fvar_scalar_template_mv,
  int_template_sv); 
  //fvar, double 
  expect_fwd_binary_scalar_matrix_eq<F, FV>(fvar_template_sv, 
  d_scalar_template_mv);  
  expect_fwd_binary_matrix_scalar_eq<F, FV>(d_scalar_template_mv, 
  fvar_template_sv);  
  expect_fwd_binary_scalar_matrix_eq<F, FV>(d_template_sv, 
  fvar_scalar_template_mv);  
  expect_fwd_binary_matrix_scalar_eq<F, FV>(fvar_scalar_template_mv, 
  d_template_sv);
  //fvar, fvar  
  expect_fwd_binary_scalar_matrix_eq<F, FV>(fvar_template_sv, 
  fvar_scalar_template_mv, 1, 0);  
  expect_fwd_binary_scalar_matrix_eq<F, FV>(fvar_template_sv, 
  fvar_scalar_template_mv, 0, 1);  
  expect_fwd_binary_scalar_matrix_eq<F, FV>(fvar_template_sv, 
  fvar_scalar_template_mv, 1, 1);  
  expect_fwd_binary_matrix_scalar_eq<F, FV>(fvar_scalar_template_mv,
  fvar_template_sv, 1, 0);  
  expect_fwd_binary_matrix_scalar_eq<F, FV>(fvar_scalar_template_mv,
  fvar_template_sv, 0, 1);  
  expect_fwd_binary_matrix_scalar_eq<F, FV>(fvar_scalar_template_mv,
  fvar_template_sv, 1, 1);  

  //vector, vector
  expect_fwd_binary_matrix_matrix_eq<F, FV>(fvar_template_mv, 
  d_template_mv);
  expect_fwd_binary_matrix_matrix_eq<F, FV>(d_template_mv, 
  fvar_template_mv);
  expect_fwd_binary_matrix_matrix_eq<F, FV>(fvar_template_mv, 
  fvar_template_mv, 1, 0);
  expect_fwd_binary_matrix_matrix_eq<F, FV>(fvar_template_mv, 
  fvar_template_mv, 0, 1);
  expect_fwd_binary_matrix_matrix_eq<F, FV>(fvar_template_mv, 
  fvar_template_mv, 1, 1);

  //scalar, vector<vector>
  //fvar, int
  expect_fwd_binary_scalar_std_vector_matrix_eq<F, FV>(int_template_sv,
  fvar_scalar_template_mv);
  expect_fwd_binary_std_vector_matrix_scalar_eq<F, FV>(
  fvar_scalar_template_mv, int_template_sv);
  //fvar, double
  expect_fwd_binary_scalar_std_vector_matrix_eq<F, FV>(fvar_template_sv,
  d_scalar_template_mv);
  expect_fwd_binary_std_vector_matrix_scalar_eq<F, FV>(d_scalar_template_mv,
  fvar_template_sv);
  expect_fwd_binary_scalar_std_vector_matrix_eq<F, FV>(d_template_sv,
  fvar_scalar_template_mv);
  expect_fwd_binary_std_vector_matrix_scalar_eq<F, FV>(
  fvar_scalar_template_mv, d_template_sv);
  //fvar, fvar
  expect_fwd_binary_scalar_std_vector_matrix_eq<F, FV>(fvar_template_sv,
  fvar_scalar_template_mv, 1, 0);
  expect_fwd_binary_scalar_std_vector_matrix_eq<F, FV>(fvar_template_sv,
  fvar_scalar_template_mv, 0, 1);
  expect_fwd_binary_scalar_std_vector_matrix_eq<F, FV>(fvar_template_sv,
  fvar_scalar_template_mv, 1, 1);
  expect_fwd_binary_std_vector_matrix_scalar_eq<F, FV>(
  fvar_scalar_template_mv, fvar_template_sv, 1, 0);
  expect_fwd_binary_std_vector_matrix_scalar_eq<F, FV>(
  fvar_scalar_template_mv, fvar_template_sv, 0, 1);
  expect_fwd_binary_std_vector_matrix_scalar_eq<F, FV>(
  fvar_scalar_template_mv, fvar_template_sv, 1, 1);

  //vector<vector>, vector<vector>
  expect_fwd_binary_std_vector_matrix_std_vector_matrix_eq<F, FV>(
  fvar_template_mv, d_template_mv);
  expect_fwd_binary_std_vector_matrix_std_vector_matrix_eq<F, FV>(
  d_template_mv, fvar_template_mv);
  expect_fwd_binary_std_vector_matrix_std_vector_matrix_eq<F, FV>(
  fvar_template_mv, fvar_template_mv, 1, 0);
  expect_fwd_binary_std_vector_matrix_std_vector_matrix_eq<F, FV>(
  fvar_template_mv, fvar_template_mv, 0, 1);
  expect_fwd_binary_std_vector_matrix_std_vector_matrix_eq<F, FV>(
  fvar_template_mv, fvar_template_mv, 1, 1);
}

#endif
