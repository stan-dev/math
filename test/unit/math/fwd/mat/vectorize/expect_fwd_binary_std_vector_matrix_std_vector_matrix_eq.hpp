#ifndef TEST_UNIT_MATH_FWD_MAT_VECTORIZE_EXPECT_FWD_BINARY_STD_VECTOR_MATRIX_STD_VECTOR_MATRIX_EQ
#define TEST_UNIT_MATH_FWD_MAT_VECTORIZE_EXPECT_FWD_BINARY_STD_VECTOR_MATRIX_STD_VECTOR_MATRIX_EQ

#include <stan/math/fwd/core/fvar.hpp>
#include <test/unit/math/fwd/mat/vectorize/build_fwd_binary_matrix.hpp>
#include <test/unit/math/fwd/mat/vectorize/expect_val_deriv_eq.hpp>
#include <vector>
#include <Eigen/Dense>

template <typename F, typename FV, typename matrix_t1, typename matrix_t2> 
void expect_fwd_binary_std_vector_matrix_std_vector_matrix_eq(
const matrix_t1& template_m1, const matrix_t2& template_m2,
bool seed_one = 1, bool seed_two = 1) {
  using std::vector;

  typedef Eigen::Matrix<FV, matrix_t1::RowsAtCompileTime, 
  matrix_t1::ColsAtCompileTime> result_matrix_t;

  const size_t num_v = 2;
  for (size_t i = 0; i < num_v; ++i) {
    for (int j = 0; j < template_m1.size(); ++j) {
      vector<matrix_t1> input_mv1;
      for (size_t k = 0; k < num_v; ++k) {
        if (seed_one && k == i)
          input_mv1.push_back(build_fwd_binary_matrix1<F>(template_m1, j));
        else
          input_mv1.push_back(build_fwd_binary_matrix1<F>(template_m1));
      }
      vector<matrix_t2> input_mv2;
      for (size_t k = 0; k < num_v; ++k) {
        if (seed_one && k == i)
          input_mv2.push_back(build_fwd_binary_matrix2<F>(template_m2, j));
        else
          input_mv2.push_back(build_fwd_binary_matrix2<F>(template_m2));
      }
      vector<result_matrix_t> fa = F::template 
      apply<vector<result_matrix_t> >(input_mv1, input_mv2);
      EXPECT_EQ(input_mv1.size(), fa.size());
      EXPECT_EQ(input_mv2.size(), fa.size());
      EXPECT_EQ(input_mv1[i].size(), fa[i].size());
      EXPECT_EQ(input_mv2[i].size(), fa[i].size());
      expect_val_deriv_eq(F::apply_base(input_mv1[i](j), input_mv2[i](j)),
      fa[i](j));
    }
  }   
}
#endif
