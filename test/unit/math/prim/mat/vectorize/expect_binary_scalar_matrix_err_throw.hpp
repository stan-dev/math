#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_BINARY_SCALAR_MATRIX_ERR_THROW_HPP
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_BINARY_SCALAR_MATRIX_ERR_THROW_HPP

#include <test/unit/math/prim/mat/vectorize/build_binary_vector.hpp>
#include <test/unit/math/prim/mat/vectorize/build_prim_binary_matrix.hpp>
#include <gtest/gtest.h>
#include <stdexcept>
#include <vector>
#include <Eigen/Dense>

template <typename F, typename result_t, typename input_t, 
typename matrix_t, int R, int C>
void expect_binary_scalar_matrix_err_throw(const std::vector<input_t>& 
input_v, const Eigen::Matrix<matrix_t, R, C>& template_m) {
  using std::vector;
  using Eigen::Matrix;

  typedef typename Eigen::Matrix<result_t, R, C> result_mt;
 
  for (size_t i = 0; i < input_v.size(); ++i) {
    matrix_t val = F::invalid_inputs2()[i];
    Eigen::Matrix<matrix_t, R, C> input_m = build_prim_binary_matrix(
    val, template_m);
    EXPECT_THROW(F::template apply<result_mt>(input_v[i], input_m), 
    std::domain_error);
  }
}

#endif
