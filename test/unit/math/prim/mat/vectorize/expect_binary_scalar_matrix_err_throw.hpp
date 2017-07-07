#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_BINARY_SCALAR_MATRIX_ERR_THROW_HPP
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_BINARY_SCALAR_MATRIX_ERR_THROW_HPP

#include <test/unit/math/prim/mat/vectorize/build_binary_vector.hpp>
#include <gtest/gtest.h>
#include <stdexcept>
#include <vector>
#include <Eigen/Dense>

template <typename F, typename result_t, typename input_t,
          typename matrix_t>
void expect_binary_scalar_matrix_err_throw(const std::vector<input_t>&
                                           input_v,
                                           const matrix_t& template_m) {
  using std::vector;
  using Eigen::Matrix;

  typedef Eigen::Matrix<result_t, matrix_t::RowsAtCompileTime,
                        matrix_t::ColsAtCompileTime> result_mt;

  for (size_t i = 0; i < input_v.size(); ++i) {
    matrix_t input_m = matrix_t::Constant(template_m.rows(),
                                          template_m.cols(),
                                          F::invalid_inputs2()[i]);
    EXPECT_THROW(F::template apply<result_mt>(input_v[i], input_m),
                 std::domain_error);
  }
}

#endif
