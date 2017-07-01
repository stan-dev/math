#ifndef TEST_UNIT_MATH_REV_MAT_VECTORIZE_EXPECT_REV_BINARY_VECTOR_VALUE_HPP
#define TEST_UNIT_MATH_REV_MAT_VECTORIZE_EXPECT_REV_BINARY_VECTOR_VALUE_HPP

#include <stan/math/rev/core/var.hpp>
#include <test/unit/math/rev/mat/vectorize/expect_rev_binary_scalar_matrix_all_eq.hpp>
#include <test/unit/math/rev/mat/vectorize/expect_rev_binary_matrix_matrix_eq.hpp>
#include <test/unit/math/rev/mat/vectorize/expect_rev_binary_scalar_std_vector_matrix_all_eq.hpp>
#include <test/unit/math/rev/mat/vectorize/expect_rev_binary_std_vector_matrix_std_vector_matrix_eq.hpp>
#include <Eigen/Dense>
#include <vector>

template <typename F, typename vector_t>
void expect_rev_binary_vector_value() {
  using stan::math::var;
  using std::vector;
  
  vector<int> int_template_sv;
  vector<double> d_template_sv;
  vector<var> var_template_sv;

  vector_t d_scalar_template_mv(5);
  vector_t d_template_mv(F::valid_inputs1().size());
  Eigen::Matrix<var, vector_t::RowsAtCompileTime, 
  vector_t::ColsAtCompileTime> var_scalar_template_mv(5);
  Eigen::Matrix<var, vector_t::RowsAtCompileTime, 
  vector_t::ColsAtCompileTime> var_template_mv(F::valid_inputs1().size());

  expect_rev_binary_scalar_matrix_all_eq<F>(int_template_sv, d_template_sv,
  var_template_sv, d_scalar_template_mv, var_scalar_template_mv);
  
  expect_rev_binary_matrix_matrix_eq<F>(var_template_mv, d_template_mv);
  expect_rev_binary_matrix_matrix_eq<F>(d_template_mv, var_template_mv);
  expect_rev_binary_matrix_matrix_eq<F>(var_template_mv, var_template_mv);

  expect_rev_binary_scalar_std_vector_matrix_all_eq<F>(
  int_template_sv, d_template_sv, var_template_sv, d_scalar_template_mv, 
  var_scalar_template_mv);

  expect_rev_binary_std_vector_matrix_std_vector_matrix_eq<F>(
  var_template_mv, d_template_mv);
  expect_rev_binary_std_vector_matrix_std_vector_matrix_eq<F>(
  d_template_mv, var_template_mv);
  expect_rev_binary_std_vector_matrix_std_vector_matrix_eq<F>(
  var_template_mv, var_template_mv);
}

#endif
