#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_BINARY_STD_VECTOR_MATRIX_SCALAR_ERR_THROW_HPP
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_BINARY_STD_VECTOR_MATRIX_SCALAR_ERR_THROW_HPP

#include <gtest/gtest.h>
#include <stdexcept>
#include <vector>
#include <Eigen/Dense>

template <typename F, typename result_t, typename input_t,
          typename matrix_t, int R, int C>
void expect_binary_std_vector_matrix_scalar_err_throw(
const Eigen::Matrix<matrix_t, R, C>& template_m,
const std::vector<input_t>& input_v) {
  using std::vector;
  using Eigen::Matrix;

  typedef typename Eigen::Matrix<matrix_t, R, C> input_mt;
  typedef typename Eigen::Matrix<result_t, R, C> result_mt;

  for (size_t i = 0; i < input_v.size(); ++i) {
    input_mt input_m = input_mt::Constant(template_m.rows(),
                                          template_m.cols(),
                                          F::invalid_inputs1()[i]);
    vector<input_mt> input_mv;
    input_mv.push_back(input_m);
    input_mv.push_back(input_m);
    EXPECT_THROW(F::template apply<vector<result_mt> >(input_mv,
                                                       input_v[i]),
                 std::domain_error);
  }
}

#endif
