#ifndef TEST_UNIT_MATH_FWD_MAT_VECTORIZE_EXPECT_FWD_BINARY_SCALAR_STD_VECTOR_MATRIX_EQ
#define TEST_UNIT_MATH_FWD_MAT_VECTORIZE_EXPECT_FWD_BINARY_SCALAR_STD_VECTOR_MATRIX_EQ

#include <stan/math/fwd/core/fvar.hpp>
#include <test/unit/math/fwd/mat/vectorize/build_binary_vector.hpp>
#include <test/unit/math/fwd/mat/vectorize/build_fwd_binary_matrix.hpp>
#include <test/unit/math/fwd/mat/vectorize/expect_val_deriv_eq.hpp>
#include <vector>
#include <Eigen/Dense>

template <typename F, typename FV, typename input_t1, typename input_t2, 
int R, int C> void expect_fwd_binary_scalar_std_vector_matrix_eq(
const std::vector<input_t1>& template_scalar_v, const Eigen::Matrix<
input_t2, R, C>& template_m, bool seed_one = 1, bool seed_two = 1) {
  using std::vector;

  typedef Eigen::Matrix<input_t2, R, C> input_matrix_t;
  typedef Eigen::Matrix<FV, R, C> result_matrix_t;

  const size_t num_v = 2;
  for (size_t i = 0; i < template_scalar_v.size(); ++i) {
    for (size_t j = 0; j < num_v; ++j) {
      for (int k = 0; k < template_m.size(); ++k) {
        input_t1 input_1; 
        if (seed_one)
          input_1 = build_binary_vector1<F>(template_scalar_v, i)[i];
        else
          input_1 = build_binary_vector1<F>(template_scalar_v)[i];
        vector<input_matrix_t> input_mv; 
        for (size_t l = 0; l < num_v; ++l) {
          if (seed_two && l == j)
            input_mv.push_back(build_fwd_binary_matrix2<F>(
            i, template_m, k));
          else
            input_mv.push_back(build_fwd_binary_matrix2<F>(i, template_m));
        }
        vector<result_matrix_t> fa = F::template 
        apply<vector<result_matrix_t> >(input_1, input_mv);
        EXPECT_EQ(input_mv.size(), fa.size());
        EXPECT_EQ(input_mv[j].size(), fa[j].size());
        expect_val_deriv_eq(F::apply_base(input_1, input_mv[j](k)), 
        fa[j](k));
      } 
    }   
  }
}
#endif
