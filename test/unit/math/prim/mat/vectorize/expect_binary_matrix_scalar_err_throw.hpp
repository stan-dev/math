#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_BINARY_MATRIX_SCALAR_ERR_THROW_HPP
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_BINARY_MATRIX_SCALAR_ERR_THROW_HPP

#include <test/unit/math/prim/mat/vectorize/build_binary_vector.hpp>
#include <gtest/gtest.h>
#include <stdexcept>
#include <vector>
#include <Eigen/Dense>

template <typename F, typename result_t, typename input_t,
          typename matrix_t>
void expect_binary_matrix_scalar_err_throw(const matrix_t& template_m, const
                                           std::vector<input_t>& input_v) {
  using std::vector;
  using Eigen::Matrix;

  typedef Eigen::Matrix<result_t, matrix_t::RowsAtCompileTime,
                        matrix_t::ColsAtCompileTime> result_mt;

  for (size_t i = 0; i < input_v.size(); ++i) {
    matrix_t input_m = matrix_t::Constant(template_m.rows(),
                                          template_m.cols(),
                                          F::invalid_inputs1()[i]);
    EXPECT_THROW(F::template apply<result_mt>(input_m, input_v[i]),
                 std::domain_error);
  }
}

#endif
