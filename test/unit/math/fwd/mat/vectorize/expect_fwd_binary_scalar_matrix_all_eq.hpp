#ifndef TEST_UNIT_MATH_FWD_MAT_VECTORIZE_EXPECT_FWD_BINARY_SCALAR_MATRIX_ALL_EQ
#define TEST_UNIT_MATH_FWD_MAT_VECTORIZE_EXPECT_FWD_BINARY_SCALAR_MATRIX_ALL_EQ

#include <stan/math/fwd/core/fvar.hpp>
#include <test/unit/math/fwd/mat/vectorize/expect_fwd_binary_scalar_matrix_eq.hpp>
#include <test/unit/math/fwd/mat/vectorize/expect_fwd_binary_matrix_scalar_eq.hpp>
#include <vector>
#include <Eigen/Dense>

template <typename F, typename FV, int R, int C> 
void expect_fwd_binary_scalar_matrix_all_eq(const std::vector<int>& 
int_template_v, const std::vector<double>& d_template_v, const
std::vector<FV> fvar_template_v, const Eigen::Matrix<double, R, C>& 
d_scalar_template_m, const Eigen::Matrix<FV, R, C>& 
fvar_scalar_template_m) {
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
#endif
